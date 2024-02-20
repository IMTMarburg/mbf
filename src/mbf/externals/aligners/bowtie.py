from .base import Aligner
import pypipegraph as ppg
import os
from pathlib import Path
import subprocess


class Bowtie(Aligner):
    @property
    def name(self):
        return "Bowtie"

    @property
    def primary_binary(self):
        return "bowtie"

    @property
    def multi_core(self):
        return True

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        return arguments + ["--threads", ncores]

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_job,
        output_bam_filename,
        parameters,
    ):
        def build_cmd():
            if hasattr(index_job, "target_folder"):  # ppg2 sharedmultifilegenjob
                index_basename = index_job.target_folder
            elif hasattr(index_job, "output_path"):  # ppg1 PrebuildJob
                index_basename = index_job.output_path
            else:
                index_basename = Path(index_job.files[0]).parent

            cmd = [
                "FROM_ALIGNER",
                self.primary_binary,
                (Path(index_basename) / "bowtie_index").absolute(),
                "-S",
            ]
            if paired_end_filename:
                cmd.extend(
                    [
                        "-1",
                        Path(paired_end_filename).absolute(),
                        "-2",
                        Path(input_fastq).absolute(),
                    ]
                )
            else:
                cmd.extend([Path(input_fastq).absolute()])
            cmd.append(str(Path(output_bam_filename).absolute()) + ".sam")
            if not "--seed" in parameters:
                parameters["--seed"] = "123123"
            for k, v in parameters.items():
                cmd.append(k)
                cmd.append(str(v))
            return cmd

        def sam_to_bam():
            import pysam

            infile = pysam.AlignmentFile(str(output_bam_filename) + ".sam", "r")
            outfile = pysam.AlignmentFile(
                str(output_bam_filename), "wb", template=infile
            )
            for s in infile:
                outfile.write(s)

        job = self.run(
            Path(output_bam_filename).parent,
            build_cmd,
            cwd=Path(output_bam_filename).parent,
            call_afterwards=sam_to_bam,
            additional_files_created=output_bam_filename,
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items()))
        )
        return job

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        if isinstance(fasta_files, (str, Path)):
            fasta_files = [fasta_files]
        if len(fasta_files) > 1:  # pragma: no cover
            raise ValueError("Bowtie can only build from a single fasta")
        cmd = [
            "FROM_ALIGNER",
            "bowtie-build",
            ",".join([str(Path(x).absolute()) for x in fasta_files]),
            (Path(output_fileprefix) / "bowtie_index").absolute(),
            "--seed",
            "123123",
        ]
        return self.get_run_func(output_fileprefix, cmd, cwd=output_fileprefix)

    def get_version(self):
        env = os.environ.copy()
        if "LD_LIBRARY_PATH" in env:  # rpy2 likes to sneak this in, breaking e.g. STAR
            del env["LD_LIBRARY_PATH"]
        res = (
            subprocess.check_output(["bowtie", "--version"], env=env)
            .decode("utf-8")
            .strip()
        )
        return res[res.find("version ") + len("version ") : res.find("\n")]

    def get_index_filenames(self):
        [
            "bowtie_index.1.ebwt",
            "bowtie_index.2.ebwt",
            "bowtie_index.3.ebwt",
            "bowtie_index.4.ebwt",
            "bowtie_index.rev.1.ebwt",
            "bowtie_index.rev.2.ebwt",
        ]


class Bowtie2(Aligner):
    @property
    def name(self):
        return "Bowtie2"

    @property
    def primary_binary(self):
        return "bowtie2"

    @property
    def multi_core(self):
        return True

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        return arguments + ["--threads", ncores]

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_job,
        output_bam_filename,
        parameters,
    ):
        def build_cmd():
            if hasattr(index_job, "target_folder"):  # ppg2 sharedmultifilegenjob
                index_basename = index_job.target_folder
            elif hasattr(index_job, "output_path"):  # ppg1 PrebuildJob
                index_basename = index_job.output_path
            else:
                index_basename = Path(index_job.files[0]).parent

            cmd = [
                "FROM_ALIGNER",
                self.primary_binary,
                "-x",
                (Path(index_basename) / "bowtie2_index").absolute(),
            ]
            if paired_end_filename:
                cmd.extend(
                    [
                        "-1",
                        Path(paired_end_filename).absolute(),
                        "-2",
                        Path(input_fastq).absolute(),
                    ]
                )
            else:
                cmd.extend([Path(input_fastq).absolute()])
            cmd.append(
                "-S",
            )
            cmd.append(str(Path(output_bam_filename).absolute()) + ".sam")
            if not "--seed" in parameters:
                parameters["--seed"] = "123123"
            for k, v in parameters.items():
                cmd.append(k)
                if v is not None:
                    cmd.append(str(v))
            return cmd

        def sam_to_bam():
            import pysam

            infile = pysam.AlignmentFile(str(output_bam_filename) + ".sam", "r")
            outfile = pysam.AlignmentFile(
                str(output_bam_filename), "wb", template=infile
            )
            for s in infile:
                outfile.write(s)

        job = self.run(
            Path(output_bam_filename).parent,
            build_cmd,
            cwd=Path(output_bam_filename).parent,
            call_afterwards=sam_to_bam,
            additional_files_created=output_bam_filename,
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items()))
        )
        return job

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        if isinstance(fasta_files, (str, Path)):
            fasta_files = [fasta_files]
        if len(fasta_files) > 1:  # pragma: no cover
            raise ValueError("Bowtie can only build from a single fasta")
        cmd = [
            "FROM_ALIGNER",
            "bowtie2-build",
            ",".join([str(Path(x).absolute()) for x in fasta_files]),
            (Path(output_fileprefix) / "bowtie2_index").absolute(),
            "--seed",
            "123123",
        ]
        return self.get_run_func(output_fileprefix, cmd, cwd=output_fileprefix)

    def get_version(self):
        res = subprocess.check_output(["bowtie2", "--version"]).decode("utf-8").strip()
        return res[res.find("version ") + len("version ") : res.find("\n")]

    def get_index_filenames(self):
        return [
            "bowtie2_index.1.bt2",
            "bowtie2_index.2.bt2",
            "bowtie2_index.3.bt2",
            "bowtie2_index.4.bt2",
            "bowtie2_index.rev.1.bt2",
            "bowtie2_index.rev.2.bt2",
        ]
