# and   import pypipegraph as ppg
#  import time
#  from pathlib import Path
#  from .. import find_code_path
from ..externals import ExternalAlgorithm


class LiftOver(ExternalAlgorithm):
    flake_name = "ucsc_tools"

    @property
    def name(self):
        return "LiftOver"

    @property
    def primary_binary(self):
        return "liftOver"

    def build_cmd(self, output_directory, ncores, arguments):  # pragma: no cover
        """Arguments = oldFile, map.chain, newFile"""
        return [self.primary_binary] + arguments

    @property
    def multi_core(self):  # pragma: no cover
        return False


class BedToBigBed(ExternalAlgorithm):
    @property
    def name(self):
        return "bedToBigBed"

    @property
    def primary_binary(self):
        return "bedToBigBed"

    def build_cmd(self, output_directory, ncores, arguments):  # pragma: no cover
        """Arguments = oldFile, map.chain, newFile"""
        return [self.primary_binary] + arguments

    @property
    def multi_core(self):  # pragma: no cover
        return False
