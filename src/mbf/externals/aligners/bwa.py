from .base import Aligner
import pypipegraph as ppg
from pathlib import Path
import subprocess
import tempfile


class BWA(Aligner):
    accepted_algos = ["mem", "aln", "samse", "sampe", "bwasw"]

    @property
    def name(self):
        return "bwa"

    @property
    def primary_binary(self):
        return "bwa"

    @property
    def multi_core(self):
        return True

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        cmd = [
            "FROM_ALIGNER",
            self.primary_binary,
            "index",
            "-p",
            str((output_fileprefix / "bwa_index").absolute()),
            "-a",
            "bwtsw",
        ]
        if not hasattr(fasta_files, "__iter__"):
            fasta_files = [fasta_files]
        cmd.extend([str(Path(x).absolute()) for x in fasta_files])
        return self.get_run_func(output_fileprefix, cmd)

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

            # output_bam_filename = Path(output_bam_filename) never used?
            algorithm = parameters.get("algorithm", "mem")
            if algorithm not in self.accepted_algos:
                raise ValueError(
                    f"Parameter 'algorithm' must be one of {self.accepted_algos}, was {algorithm}."
                )
            cmd = ["FROM_ALIGNER", str(self.path / "bwa-0.7.17" / "bwa"), algorithm]
            if algorithm == "mem":
                cmd.extend(["-k", str(parameters.get("-k", 19))])
                cmd.extend(["-w", str(parameters.get("-w", 100))])
                cmd.extend(["-d", str(parameters.get("-d", 100))])
                cmd.extend(["-r", str(parameters.get("-r", 1.5))])
                cmd.extend(["-c", str(parameters.get("-c", 10000))])
                cmd.extend(["-A", str(parameters.get("-A", 1))])
                cmd.extend(["-B", str(parameters.get("-B", 4))])
                cmd.extend(["-O", str(parameters.get("-O", 6))])
                cmd.extend(["-E", str(parameters.get("-E", 1))])
                cmd.extend(["-L", str(parameters.get("-L", 5))])
                cmd.extend(["-U", str(parameters.get("-U", 9))])
                cmd.extend(["-T", str(parameters.get("-T", 30))])
            cmd.append(str(index_basename / "bwa_index"))
            cmd.append(str(Path(input_fastq).absolute()))
            if paired_end_filename:
                cmd.append(str(Path(paired_end_filename).absolute()))
            return cmd

        job = self.run(
            Path(output_bam_filename).parent,
            build_cmd,
            additional_files_created=[
                output_bam_filename,
            ],
            call_afterwards=self.sam_to_bam(
                output_bam_filename.parent / "stdout.txt", output_bam_filename
            ),
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items()))
        )
        return job

    def sam_to_bam(self, infile: Path, outfile: Path):
        def __convert():
            cmd = ["samtools", "view", "-b", str(infile)]
            tmp = tempfile.NamedTemporaryFile("w")
            with tmp:
                subprocess.check_call(cmd, stdout=tmp)
                cmd = ["samtools", "sort", tmp.name]
                subprocess.check_call(cmd, stdout=outfile.open("w"))
            with infile.open("w") as op:
                op.write("Moved to bam via samtools")

        return __convert

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        if "mem" == arguments[1]:
            return arguments + ["-t", str(ncores)]
        else:
            return arguments

    def get_index_filenames(self):
        return [
            "bwa_index.amb",
            "bwa_index.ann",
            "bwa_index.bwt",
            "bwa_index.pac",
            "bwa_index.sa",
        ]
