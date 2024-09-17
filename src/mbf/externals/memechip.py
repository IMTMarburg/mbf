from .externals import ExternalAlgorithm
import shutil
from pathlib import Path
import pypipegraph as ppg


class MemeChIP(ExternalAlgorithm):
    @property
    def name(self):
        return "memechip"

    def run_on_gr(self, gr, options=None):
        output_directory = gr.result_dir / "MemeChIP"
        if output_directory.exists():
            shutil.rmtree(output_directory)
        # memechip does this output_directory.mkdir(exist_ok=True, parents=True)

        input = gr.to_fasta()
        args = [
            "-oc",
            str(output_directory.absolute()),
        ]
        if options:
            args.extend(options)

        args.append(str(input.files[0].absolute()))

        res = self.run(
            output_directory,
            # fmt: off
            args,
            # fmt: on
            additional_files_created=[
                output_directory / "meme-chip.html",
                output_directory / "summary.tsv",
                output_directory / "combined.meme",
            ],
        )
        res.depends_on(input)
        return res

    def build_cmd(self, output_directory, ncores, arguments):  # pragma: no cover
        return [self.primary_binary] + arguments

    @property
    def multi_core(self):  # pragma: no cover
        return False

    @property
    def primary_binary(self):  # pragma: no cover
        return "meme-chip"
