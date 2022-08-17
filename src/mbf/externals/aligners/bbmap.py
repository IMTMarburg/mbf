from .base import Aligner

try:
    import mbf.align as mbf_align
except ImportError:
    mbf_align = None
from pathlib import Path
import subprocess


class BBMap(Aligner):
    @property
    def name(self):
        return "BBMap"

    @property
    def primary_binary(self):
        return "bbmap.sh"

    @property
    def multi_core(self):
        return False

    # todo: flake me
    # latest version "38.86"
    # def fetch_version(self, version, target_filename):  # pragma: no cover
    # url = f"https://sourceforge.net/projects/bbmap/files/BBMap_{version}.tar.gz/download"
    #        cmd = ["curl", url, "-L", "--output", "bb.tar.gz"]
    #        subprocess.check_call(cmd)
    # with open(target_filename, "wb") as op:
    # download_file(url, op)

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        raise NotImplementedError

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_job,
        output_bam_filename,
        parameters,
    ):
        raise NotImplementedError

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        raise NotImplementedError

    def get_index_filenames(self):
        raise NotImplementedError("Find outh what files bbmap leaves")


if mbf_align:

    class ExtendCigarBBMap(mbf_align.post_process._PostProcessor):
        def __init__(self, samformat="1.4"):
            self.samformat = samformat
            self.bbmap = BBMap()
            self.name = "BBMap_reformat"
            self.result_folder_name = Path("results") / "aligned" / "ExtendCigarBBMap"

        def process(self, input_bam_name, output_bam_name, result_dir):
            cmd = [
                str(
                    self.bbmap.path / "bbmap" / "reformat.sh"
                ),  # todo: flake / nix variant
                f"in={str(input_bam_name.absolute().resolve())}",
                f"out={str(output_bam_name.absolute().resolve())}",
                f"sam={self.samformat}",
            ]
            print(" ".join(cmd))
            subprocess.check_call(cmd)

        def register_qc(self, new_lane):
            pass  # pragma: no cover

        def get_version(self):
            return self.bbmap.version

        def get_parameters(self):
            return self.samformat
