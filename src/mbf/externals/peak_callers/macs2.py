"""
This is a wrapper for MACS2 analysis
"""

from ..externals import ExternalAlgorithm
from pathlib import Path
import pypipegraph as ppg


class MACS2(ExternalAlgorithm):
    """Peak calling for model based analysis of chip seq data"""

    @property
    def name(self):
        return "MACS2"

    @property
    def primary_binary(self):
        return "macs2"

    @property
    def multi_core(self):
        return False

    def build_cmd(self, output_directory, ncores, arguments):
        cmd = [
            self.primary_binary,
            "callpeak",
        ]
        cmd.append("-t")
        cmd.extend(arguments["input_bam"])
        if "background_bam" in arguments:
            cmd.append("-c")
            cmd.extend(arguments["background_bam"])
        cmd.extend(["-f", "BAM"])
        cmd.extend(["--outdir", output_directory])
        for param, value in arguments["parameters"].items():
            if value is None or value == "":
                cmd.append(param)
            else:
                cmd.append(f"{param}={value}")

        return cmd

    def call_peaks(
        self,
        input_lane,
        background_lane=None,
        parameters={},
        name=None,
        result_dir=None,
    ):
        """Call peaks using macs.
        Use {"--nomodel": False} for no-value parameters
        """
        from mbf.genomics.regions import GenomicRegions  # avoid circular imports?

        if not isinstance(input_lane, list):
            input_lane = [input_lane]

        if background_lane is not None:
            if not isinstance(background_lane, list):
                background_lane = [background_lane]
            if name is None:
                name = f"{input_lane[0].name}_vs_{background_lane[0].name}_MACS2"
            args = {
                "input_bam": [x.get_bam_names()[0] for x in input_lane],
                # "paired": input_lane.is_paired,
                "background_bam": [x.get_bam_names()[0] for x in background_lane],
                "parameters": parameters,
            }
        else:
            if name is None:
                name = f"{input_lane[0].name}_MACS2_{self.version}"
            args = {
                "input_bam": (x.get_bam_names()[0] for x in input_lane),
                # "paired": input_lane.is_paired,
                "parameters": parameters,
            }

        output_directory = Path("results", "MACS2", name)
        add_files = [output_directory / "NA_peaks.xls"]
        if '--bdg' in parameters:
            add_files.append(output_directory / "NA_treat_pileup.bdg")
            add_files.append(output_directory / "NA_control_lambda.bdg")
        run_job = self.run(
            output_directory,
            args,
            additional_files_created=add_files,
        )
        run_job.depends_on(
            ppg.ParameterInvariant(name, parameters),
            [x.load() for x in input_lane],
        )
        if background_lane is not None:
            run_job.depends_on(
                [x.load() for x in background_lane],
            )

        def do_load_peaks():
            import pandas as pd

            df = pd.read_csv(
                output_directory / "NA_peaks.xls", sep="\t", comment="#", skiprows=16
            )
            df["chr"] = [str(x) for x in df["chr"]]
            df["start"] -= 1  # correct for macs one base offset
            df["start"] = df["start"].clip(lower=0)
            df["end"] -= 1
            renames = {"end": "stop"}
            for f in [
                "-log10(pvalue)",
                "fold_enrichment",
                "-log10(qvalue)",
            ]:
                renames[f] = "MACS2 " + f
            df = df.rename(columns=renames)
            fdr = [pow(10, -x) for x in df["MACS2 -log10(qvalue)"]]
            df["FDR"] = fdr
            df = df.drop("name", axis=1)
            return df

        if background_lane is not None:
            vid = (
                [x.vid for x in input_lane] + ["vs"] + [x.vid for x in background_lane]
            )
        else:
            vid = [x.vid for x in input_lane]

        return GenomicRegions(
            name,
            do_load_peaks,
            run_job,
            input_lane[0].genome,
            sheet_name="Peaks",
            # vid=[input_lane.vid, "vs", background_lane.vid],
            on_overlap="ignore",
            result_dir=result_dir,
            vid=vid,
        )

    def bdgcmp(self, treatment_peaks, parameters=[], name=None, result_dir=None):
        if name is None:
            raise ValueError("Name is required")
            # this never worked...
            # name = f"{treatment_lane.name}_against_{control_lane.name}"

        output_directory = Path("results", "MACS2", "bdgcmp", name)
        output_directory.mkdir(parents=True, exist_ok=True)

        treatment_file = f"results/MACS2/{treatment_peaks.name}/NA_treat_pileup.bdg"
        control_file = f"results/MACS2/{treatment_peaks.name}/NA_control_lambda.bdg"

        cmd = [
            "bdgcmp",
        ]
        cmd.extend(["-t", treatment_file])
        cmd.extend(["-c", control_file])
        cmd.extend(["--outdir", output_directory])
        cmd.extend(["--o-prefix", name])
        for param in parameters:
            cmd.append(param)

        run_job = self.run(
            output_directory,
            cmd,
        )

        run_job.depends_on(
            ppg.ParameterInvariant(name, parameters),
            treatment_peaks.load(),
        )
        run_job.depends_on_file(treatment_file)
        run_job.depends_on_file(control_file)

        return run_job
