from ..externals import ExternalAlgorithm
from pathlib import Path
import pypipegraph as ppg


class PeakZilla(ExternalAlgorithm):
    @property
    def name(self):
        return "PeakZilla"

    @property
    def primary_binary(self):
        return "peakzilla.py"

    @property
    def multi_core(self):
        return False

    def call_peaks(self, input_lane, background_lane, parameters={}, name=None):
        from mbf.genomics.regions import GenomicRegions

        if input_lane.genome != background_lane.genome:  # pragma: no cover
            raise ValueError("Lanes had unequal genomes")

        if name is None:
            name = (
                f"{input_lane.name}_vs_{background_lane.name}_Peakzilla_{self.version}"
            )

        output_directory = Path("cache", "PeakZilla", name)
        run_job = self.run(
            output_directory,
            {
                "input_bam": input_lane.get_bam_names()[0],
                "paired": input_lane.is_paired,
                "background_bam": background_lane.get_bam_names()[0],
                "parameters": parameters,
            },
        )
        run_job.depends_on(
            ppg.ParameterInvariant(name, parameters),
            input_lane.load(),
            background_lane.load(),
        )

        def do_load_peaks():
            import pandas as pd

            df = pd.read_csv(output_directory / "stdout.txt", sep="\t")
            df = df.rename(
                columns={
                    "#Chromosome": "chr",
                    "Start": "start",
                    "End": "stop",
                    "Name": "pz-Name",
                    "Summit": "pz-Summit",
                    "Score": "pz-Score",
                    "ChIP": "pz-chIP",
                    "Control": "pz-Control",
                    "FoldEnrichment": "pz-FoldEnrichment",
                    "DistributionScore": "pz-DistributionScore",
                    "FDR": "pz-FDR",
                }
            )
            df["chr"] = df["chr"].astype(str)
            df["start"] -= 1
            df["start"] = df["start"].clip(lower=0)
            df["stop"] -= 1
            return df

        return GenomicRegions(
            name,
            do_load_peaks,
            run_job,
            input_lane.genome,
            sheet_name="Peaks",
            vid=[input_lane.vid, "vs", background_lane.vid],
            on_overlap="ignore",
        )

    def build_cmd(self, output_directory, ncores, arguments):
        cmd = [
            "peakzilla.py",
            "-l",
            str(output_directory / "log.txt"),
        ]
        if arguments["paired"]:  # pragma: no cover
            cmd.append("-p")
        cmd.append(arguments["input_bam"])
        cmd.append(arguments["background_bam"])
        for k, v in arguments["parameters"].items():
            cmd.append(k)
            cmd.append(v)
        return cmd
