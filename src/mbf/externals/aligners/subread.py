from pypipegraph.util import output_file_exists
from .base import Aligner
import pypipegraph as ppg
from pathlib import Path
from ..util import Version
import subprocess


class Subread(Aligner):
    @property
    def name(self):
        return "Subread"

    @property
    def primary_binary(self):
        return "subread-align"

    @property
    def multi_core(self):
        return True

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        if "subread-align" in arguments[0]:
            return arguments + ["-T", str(ncores)]
        else:
            return arguments

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_job,
        output_bam_filename,
        parameters,
    ):
        output_bam_filename = Path(output_bam_filename)
        if not parameters.get("input_type") in ("dna", "rna"):
            raise ValueError("invalid parameters['input_type'], must be dna or rna")
        if isinstance(index_job, Path):
            index_job = [
                ppg.FileInvariant(index_basename / x)
                for x in self.get_index_filenames()
            ]

        def build_cmd():
            if hasattr(index_job, "target_folder"):  # ppg2 sharedmultifilegenjob
                index_basename = index_job.target_folder
            elif hasattr(index_job, "output_path"):  # ppg1 PrebuildJob
                index_basename = index_job.output_path
            else:
                index_basename = Path(index_job.files[0]).parent

            if parameters["input_type"] == "dna":
                input_type = "1"
            else:
                input_type = "0"
            p_output_bam_filename = Path(output_bam_filename)
            cmd = [
                "FROM_ALIGNER",
                "subread-align",
                "-t",
                input_type,
                "-I",
                "%i" % parameters.get("indels_up_to", 5),
                "-M",
                "%i" % parameters.get("-M", 3),
                "-B",
                "%i" % parameters.get("max_mapping_locations", 1),
                "-i",
                (Path(index_basename) / "subread_index").absolute(),
                "-r",
                Path(input_fastq).absolute(),
                "-o",
                p_output_bam_filename.absolute(),
            ]
            if "keepReadOrder" in parameters:  # pragma: no cover
                cmd.append("--keepReadOrder")
            else:
                cmd.append("--sortReadsByCoordinates")
            if paired_end_filename:
                cmd.extend(("-R", str(Path(paired_end_filename).absolute())))
            if "max_mapping_locations" in parameters:  # pragma: no cover
                cmd.append("--multiMapping")
            if "--sv" in parameters:
                cmd.append("--sv")
            if "-m" in parameters:
                cmd.extend(["-m", parameters["-m"]])
            return cmd

        def remove_bai():
            # subread create broken bais where idxstat doesn't work.
            # but the mbf_aligned.lanes.AlignedSample will recreate it.
            # so it's ok if we simply throw it away here.
            bai_name = output_bam_filename.with_name(output_bam_filename.name + ".bai")
            if bai_name.exists():
                bai_name.unlink()

        job = self.run(
            Path(output_bam_filename).parent,
            build_cmd,
            additional_files_created=[
                output_bam_filename,
                output_bam_filename.with_name(output_bam_filename.name + ".indel.vcf"),
                # subread create broken bais where idxstat doesn't work.
                # output_bam_filename.with_name(output_bam_filename.name + ".bai"),
            ],
            call_afterwards=remove_bai,
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items()))
        )
        job.depends_on(index_job)
        return job

    def get_version(self):
        p = subprocess.Popen(
            ["subread-align", "-version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = p.communicate()
        res = stderr.decode("utf-8").strip()
        if not "Subread-align v" in res:
            raise ValueError("Failed to parse version", stderr)
        return res[res.rfind("v") + 1 :]

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        cmd = [
            "FROM_ALIGNER",
            "subread-buildindex",
            "-o",
            str((output_fileprefix / "subread_index").absolute()),
        ]
        if not hasattr(fasta_files, "__iter__"):
            fasta_files = [fasta_files]
        cmd.extend([str(Path(x).absolute()) for x in fasta_files])
        return self.get_run_func(output_fileprefix, cmd)

    def get_index_version_range(self):
        """What minimum_acceptable_version, maximum_acceptable_version for the index is ok?"""
        if Version(self.version) >= "1.6":
            return "1.6", None
        else:
            return "0.1", "1.5.99"

    def get_index_filenames(self):
        return [
            "subread_index.00.b.array",
            "subread_index.00.b.tab",
            "subread_index.files",
            "subread_index.reads",
        ]

    def get_alignment_stats(self, output_bam_filename):
        import re

        output_bam_filename = Path(output_bam_filename)
        target = output_bam_filename.parent / "stderr.txt"
        raw = target.read_text()
        raw = raw[raw.find("Summary =") :]
        result = {}
        keys = "Uniquely mapped", "Multi-mapping", "Unmapped"
        for k in keys:
            try:
                result[k] = int(re.findall(f"{k} : (\\d+)", raw)[0])
            except IndexError as e:
                raise KeyError(k, e)
        return result


class Subjunc(Subread):
    @property
    def name(self):
        return "Subjunc"

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        if "subjunc" in arguments[0]:
            return arguments + ["-T", str(ncores)]
        else:
            return arguments

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_basename,
        output_bam_filename,
        parameters,
    ):
        output_bam_filename = Path(output_bam_filename)
        cmd = [
            "FROM_ALIGNER",
            "subjunc",
            "-I",
            "%i" % parameters.get("indels_up_to", 5),
            "-M",
            "%i" % parameters.get("-M", 3),
            "-B",
            "%i" % parameters.get("max_mapping_locations", 1),
            "-i",
            (Path(index_basename) / "subread_index").absolute(),
            "-r",
            Path(input_fastq).absolute(),
            "-o",
            output_bam_filename.absolute(),
        ]
        if "keepReadOrder" in parameters:  # pragma: no cover
            cmd.append("--keepReadOrder")
        else:
            cmd.append("--sortReadsByCoordinates")
        if paired_end_filename:
            cmd.extend(("-R", str(Path(paired_end_filename).absolute())))
        if "max_mapping_locations" in parameters:  # pragma: no cover
            cmd.append("--multiMapping")
        if "-m" in parameters:
            cmd.extend(["-m", parameters["-m"]])
        if "--complexIndels" in parameters:
            cmd.append("--complexIndels")

        def remove_bai():
            # subread create broken bais where idxstat doesn't work.
            # but the mbf_aligned.lanes.AlignedSample will recreate it.
            # so it's ok if we simply throw it away here.
            bai_name = output_bam_filename.with_name(output_bam_filename.name + ".bai")
            if bai_name.exists():
                bai_name.unlink()

        job = self.run(
            Path(output_bam_filename).parent,
            cmd,
            additional_files_created=[
                output_bam_filename,
                # subread create broken bais where idxstat doesn't work.
                # output_bam_filename.with_name(output_bam_filename.name + ".bai"),
            ],
            call_afterwards=remove_bai,
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items()))
        )
        return job
