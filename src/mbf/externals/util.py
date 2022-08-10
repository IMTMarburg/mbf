import functools
import tempfile
import subprocess
import os
from pathlib import Path
from contextlib import contextmanager


class lazy_property(object):
    """
    meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    """

    def __init__(self, fget):
        self.fget = fget

        # copy the getter function's docstring and other attributes
        functools.update_wrapper(self, fget)

    def __get__(self, obj, cls):
        # if obj is None: # this was in the original recepie, but I don't see
        # when it would be called?
        # return self

        value = self.fget(obj)
        setattr(obj, self.fget.__name__, value)
        return value


def lazy_method(func):
    """
    meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    """
    cache_name = "_cached_" + func.__name__

    def inner(self):
        if not hasattr(self, cache_name):
            setattr(self, cache_name, func(self))
        return getattr(self, cache_name)

    return inner


def sort_versions(versions):
    """Sort versions, from natsort manual:
    Sorts like this:
        ['1.1', '1.2', '1.2alpha', '1.2beta1', '1.2beta2', '1.2rc1', '1.2.1', '1.3']
    """
    import natsort

    return natsort.natsorted(
        versions,
        key=lambda x: x.replace(".", "~")
        if not isinstance(x, tuple)
        else x[0].replace(".", "~"),
    )


class Version:
    """A smartish comparable version tuple.

    used e.g. do decide whether an aligners index is
    ok to use for the given aligner version

    """

    def __init__(self, version):
        self.version = version

    def __str__(self):
        return self.version

    def __repr__(self):
        return 'Version("%s")' % (self.version,)

    def __eq__(self, other_version):
        if isinstance(other_version, Version):
            other = other_version.version
        else:
            other = other_version
        return self.version == other

    def __lt__(self, other_version):
        if isinstance(other_version, Version):
            other = other_version.version
        else:
            other = other_version
        s = sort_versions([self.version, other])
        return s[0] == self.version and not self.version == other

    def __le__(self, other_version):
        return (self == other_version) or (self < other_version)

    def __gt__(self, other_version):
        if isinstance(other_version, Version):
            other = other_version.version
        else:
            other = other_version
        s = sort_versions([self.version, other])
        return s[1] == self.version and not self.version == other

    def __ge__(self, other_version):
        return (self == other_version) or (self > other_version)


# ppg1 only
class UpstreamChangedError(ValueError):
    pass


def download_file(url, file_object):
    """Download an url"""
    if isinstance(file_object, (str, Path)):
        raise ValueError("download_file needs a file-object not a name")

    try:
        if url.startswith("ftp"):
            return download_ftp(url, file_object)
        else:
            return download_http(url, file_object)
    except Exception as e:
        raise ValueError("Could not download %s, exception: %s" % (repr(url), e))


def download_http(url, file_object):
    """Download a file from http"""
    import requests
    import shutil

    r = requests.get(url, stream=True)
    if r.status_code != 200:
        raise ValueError("HTTP Error return: %i fetching %s" % (r.status_code, url))
    r.raw.decode_content = True
    shutil.copyfileobj(r.raw, file_object)


def download_ftp(url, file_object):
    """Download a file from ftp"""
    import ftplib
    import urllib

    schema, host, path, parameters, query, fragment = urllib.parse.urlparse(url)
    if "ftp_proxy" in os.environ:
        tf = tempfile.NamedTemporaryFile()
        cmd = ["wget", "-O", tf.name, url]
        subprocess.check_call(cmd)
        tf.seek(0, 0)
        bs = 1024 * 1024 * 10
        block = tf.read(bs)
        while block:
            file_object.write(block)
            block = tf.read(bs)
        tf.close()

    else:
        with ftplib.FTP(host) as ftp:
            try:
                ftp.set_pasv(True)
                ftp.login("anonymous", "")
                if "\n" in path:  # pragma: no cover
                    raise ValueError("New line in path: %s" % (repr(path),))
                if path.endswith("/"):
                    ftp.retrbinary("LIST " + path, file_object.write)
                else:
                    ftp.retrbinary("RETR " + path, file_object.write)
            except ftplib.Error as e:
                raise ValueError("Error retrieving urls %s: %s" % (url, e))


def get_page(url):
    """Download a web page (http/ftp) into a string"""
    from io import BytesIO

    tf = BytesIO()
    download_file(url, tf)
    tf.seek(0, 0)
    return tf.read().decode("utf-8")


def download_file_and_gunzip(url, unzipped_filename):
    import shutil
    import gzip
    import tempfile

    tf = tempfile.NamedTemporaryFile(suffix=".gz")
    download_file(url, tf)
    tf.flush()

    with gzip.GzipFile(tf.name, "rb") as gz_in:
        with open(unzipped_filename, "wb") as op:
            shutil.copyfileobj(gz_in, op)


def download_file_and_gzip(url, gzipped_filename):
    import shutil
    import gzip
    import tempfile

    gzipped_filename = str(gzipped_filename)
    if not gzipped_filename.endswith(".gz"):  # pragma: no cover
        raise ValueError("output filename did not end with .gz")

    with tempfile.NamedTemporaryFile(suffix="") as tf:
        with gzip.GzipFile(tf.name, "wb") as gf:
            download_file(url, gf)
        shutil.copy(tf.name, gzipped_filename)


def write_md5_sum(filepath):
    """Create filepath.md5sum with the md5 hexdigest"""
    from pypipegraph.util import checksum_file

    md5sum = checksum_file(filepath)
    (filepath.with_name(filepath.name + ".md5sum")).write_text(md5sum)


def to_string(s, encoding="utf-8"):
    if isinstance(s, str):
        return s
    else:
        return s.decode(encoding)


def to_bytes(x, encoding="utf-8"):
    """In python3: str -> bytes. Bytes stay bytes"""
    if isinstance(x, bytes):
        return x
    else:
        return x.encode(encoding)


def chmod(filename, mode):
    """Chmod if possible - otherwise try to steal the file and chmod then"""
    import os
    import shutil

    try:
        os.chmod(filename, mode)
    except OSError as e:  # pragma: no cover
        if (
            str(e).find("Operation not permitted") == -1
            and str(e).find("Permission denied") == -1
        ):
            raise
        else:  # steal ownership and set the permissions...
            t = filename + ".temp"
            shutil.copyfile(filename, t)
            try:
                os.chmod(t, mode)
            except OSError:
                pass
            os.unlink(filename)
            shutil.move(t, filename)


def download_zip_and_turn_into_tar_gzip(url, target_filename, chmod_x_files=[]):
    """Download a zip archive and turn it into the correct tar.gzip"""
    import tempfile
    import subprocess
    from .externals import reproducible_tar

    if isinstance(chmod_x_files, str):  # pragma: no cover
        chmod_x_files = [chmod_x_files]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        with (tmpdir / "source.zip").open("wb") as zip_file:
            download_file(url, zip_file)
        import zipfile

        with zipfile.ZipFile(zip_file.name, "r") as zip_ref:
            zip_ref.extractall(tmpdir / "target")
        for fn in chmod_x_files:
            subprocess.check_call(
                ["chmod", "+x", str(tmpdir.absolute() / "target" / fn)]
            )
        reproducible_tar(target_filename.absolute(), "./", cwd=tmpdir / "target")


def download_mercurial_update_and_zip(url, changeset, target_filename):
    """Download a mercurial repo, update it to a specific changeset, and tar it up"""
    import tempfile
    import subprocess
    from .externals import reproducible_tar

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        subprocess.check_call(["hg", "clone", url, str(tmpdir.absolute())])
        subprocess.check_call(["hg", "up", "-r", changeset], cwd=tmpdir)
        reproducible_tar(target_filename.absolute(), "./", cwd=tmpdir)


def download_tar_bz2_and_turn_into_tar_gzip(
    url, target_filename, version, chmod_x_files=[], make=True
):
    """Download a tar.bz2 archive and turn it into the correct tar.gzip"""

    import tempfile
    import subprocess
    from .externals import reproducible_tar
    import tarfile

    if isinstance(chmod_x_files, str):  # pragma: no cover
        chmod_x_files = [chmod_x_files]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        with (tmpdir / "source.tar.bz2").open("wb") as archive_file:
            download_file(url, archive_file)

        with tarfile.open(archive_file.name, "r") as zip_ref:
            zip_ref.extractall(tmpdir / "target")
        for fn in chmod_x_files:
            subprocess.check_call(
                ["chmod", "+x", str(tmpdir.absolute() / "target" / fn)]
            )
        if make:

            wd = [x for x in (tmpdir.absolute() / "target").iterdir()][0]
            subprocess.check_call(["make"], cwd=wd)
        reproducible_tar(target_filename.absolute(), "./", cwd=tmpdir / "target")


def binary_exists(binary_name):
    import subprocess
    import shlex

    p = subprocess.Popen(
        f"command -v {shlex.quote(binary_name)}", shell=True, stdout=subprocess.PIPE
    )
    p.communicate()
    return p.returncode == 0



@contextmanager
def chdir(path: Path):
    """change the current working directory, 
    in a context manager
    """

    original = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(original)
