import bz2
import subprocess
import tarfile
import gzip
import os


def read_fastq_iterator(file_object):
    """A very dump and simple fastq reader, mostly for testing the other more sophisticated variants

    Yield (seq, name, quality)
    """
    row1 = file_object.readline()
    row2 = file_object.readline()
    row3 = file_object.readline()
    row4 = file_object.readline()
    while row1:
        seq = row2[:-1]
        quality = row4[:-1]
        name = row1[1:-1]
        yield (seq, name, quality)
        row1 = file_object.readline()
        row2 = file_object.readline()
        row3 = file_object.readline()
        row4 = file_object.readline()


class BlockedFileAdaptor:
    """Wraps read_file_blocked into a semi file like object, good enough for read_fastq_iterator"""

    def __init__(self, filename, block_size=1 * 1024 * 1024):
        self.iterator = None
        self.name = filename
        self.line_iterator = None
        self.block_size = block_size

    def readline(self):
        if self.line_iterator is None:
            self.line_iterator = self._iter_lines()
        return next(self.line_iterator)

    def _iter_lines(self):
        if not self.iterator:
            self.iterator = read_file_blocked(self.name, self.block_size)
        remainder = b""
        for block in self.iterator:
            block = remainder + block
            lines = block.split(b"\n")
            if block.endswith(b"\n"):
                remainder = b""
                del lines[-1]  # remove the empty line...
            else:
                remainder = lines[-1]
                del lines[-1]
            for l in lines:
                yield l + b"\n"
        if remainder:
            yield remainder


def read_file_blocked(filename, block_size=1 * 1024 * 1024):
    """Read a (possibly compressed) file block by uncompressed block, yielding the blocks"""
    filename = str(filename)  # Work with Paths...
    p = None
    zf = None
    if not os.path.exists(filename):
        raise IOError("File not found: %s" % filename)
    if filename.endswith(".tar.bz2"):
        cmd = ["lbzip2 -d | tar -xO"]
        stdin = open(filename, "rb")
        p = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=stdin, shell=True
        )
    elif filename.endswith(".gz"):
        cmd = ["gzip", "-cd", filename]
        p = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False
        )
    elif filename.endswith(".tgz") or filename.endswith(".tar.gz"):
        if not os.path.exists(filename):
            raise IOError("[Errno 2] No such file or directory: '%s" % filename)
        cmd = [
            "tar",
            "-xOf",
            filename,
        ]  # as cool as python is, 5x speedup by calling tar directly instead of using tarfile is *not* funny.
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            bufsize=0,
        )
    elif filename.endswith(".bz2"):
        op = bz2.BZ2File(filename)
    elif filename.endswith(".zip"):
        import zipfile

        zf = zipfile.ZipFile(filename, "r")
        names = zf.namelist()
        if len(names) > 1:
            raise ValueError(
                "This zipfile contains more than one compressed sequence file. That's unexpected, augment read_file_blocked"
            )
        op = zf.open(names[0])
    elif filename.endswith(".tar"):
        # Note delegation instead of using the common control flow... since it generates multiple pipes
        for block in _read_file_blocked_tarfile(filename, block_size):
            yield block
        return  # no need for all the rest...

    else:
        op = open(filename, "rb")
    if not p is None:
        input_pipe = p.stdout
        if not p.stdin is None:
            p.stdin.close()
    else:
        input_pipe = op
    try:
        while True:
            block = input_pipe.read(block_size)
            if not block:
                break
            else:
                yield block
                del block
    finally:
        input_pipe.close()
        if not p is None:
            p.wait()
            p.stderr.close()
        if not zf is None:
            zf.close()


def _read_file_blocked_tarfile(filename, block_size):
    # assume for now that it's a tar of multiple fastq.gz, and possibly a csv that's being ignored, and will explode on other files... Such data is for example produced by Braunschweig after their 1.8 Illumina upgrade
    try:
        tf = tarfile.open(filename)
        for tarinfo in tf.getmembers():
            if tarinfo.size == 0:  # ignore empty entries (directories)
                continue
            if tarinfo.name.endswith(".csv"):
                continue
            if tarinfo.name.endswith(".fastq.gz"):
                cmd = "tar -xvOf %s %s | gzip -cd" % (filename, tarinfo.name)
                p = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    bufsize=0,
                    shell=True,
                )
                try:
                    while True:
                        block = p.stdout.read(block_size)
                        if not block:
                            break
                        else:
                            yield block
                            del block
                finally:
                    p.wait()
                    p.stdout.close()
                    p.stderr.close()
            else:
                raise ValueError(
                    "chipseq.common.read_file_blocked does not know how to handle this file (%s) in tar archive %s"
                    % (tarinfo.name, filename)
                )
    finally:
        tf.close()
