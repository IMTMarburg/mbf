# and   import pypipegraph as ppg
#  import time
#  from pathlib import Path
#  from .. import find_code_path
from .externals import ExternalAlgorithm
from pathlib import Path
import pypipegraph as ppg


class Snpdiffrs(ExternalAlgorithm):
    @property
    def name(self):
        return "snpdiffrs"

    def run_n_to_n(
        self,
        output_directory,
        sample_to_bams: dict,
        min_score=None,
        quality_threshold=None,
        filter_homo_polymer_threshold=None,
        chromosomes=None,
    ):
        output_directory = Path(output_directory)
        output_directory.mkdir(exist_ok=True, parents=True)
        for name, bams in sample_to_bams.items():
            if not isinstance(bams, list):
                raise ValueError("Sample-bams must be a list per sample")

        def write_toml(output_filename):
            toml = [f"output_dir = '{str(output_directory)}'"]
            if min_score:
                toml.append(f"min_score = {min_score:.2f}")
            if quality_threshold:
                toml.append(f"quality_threshold = {int(quality_threshold)}")
            if filter_homo_polymer_threshold:
                toml.append(
                    f"filter_homo_polymer_threshold = {int(filter_homo_polymer_threshold)}"
                )
            if chromosomes:
                toml.append(f"chromosomes = {chromosomes}")
            toml.append("[samples]")
            for sample_name, bams in sample_to_bams.items():
                toml.append(f"\t{sample_name} = {[str(x) for x in bams]}")
            Path(output_filename).write_text("\n".join(toml))

        prep_job = ppg.FileGeneratingJob(
            output_directory / "input.toml", write_toml
        ).depends_on(
            ppg.ParameterInvariant(
                output_directory / "input.toml",
                (
                    sample_to_bams,
                    min_score,
                    quality_threshold,
                    filter_homo_polymer_threshold,
                    chromosomes,
                ),
            )
        )
        res = self.run(output_directory, [str(output_directory / "input.toml")])
        res.depends_on(prep_job)
        return res

    def build_cmd(self, output_directory, ncores, arguments):  # pragma: no cover
        return [self.primary_binary] + arguments

    @property
    def multi_core(self):  # pragma: no cover
        return True

    @property
    def primary_binary(self):  # pragma: no cover
        return "snpdiffrs"
