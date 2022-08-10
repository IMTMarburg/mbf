#  import pypipegraph as ppg
#  import time
#  from pathlib import Path
#  from .. import find_code_path
from ..externals import ExternalAlgorithm


class FASTQC(ExternalAlgorithm):
    @property
    def primary_binary(self):
        return "fastqc"

    @property
    def name(self):
        return "FASTQC"

    @property
    def multi_core(self):
        return False  # fastqc has a threads option - and does not make use of it

    def build_cmd(self, output_directory, ncores, arguments):
        input_files = arguments
        return [
            self.primary_binary,
            "-t",
            str(ncores),
            "--noextract",
            "--quiet",
            "-o",
            str(output_directory),
        ] + [str(x) for x in input_files]
