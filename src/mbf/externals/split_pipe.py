from pkg_resources import add_activation_listener
import pypipegraph2 as ppg
import os
import subprocess
from pathlib import Path

from ..externals import ExternalAlgorithm

_genome_job_cache = {}


class SplitPipe(ExternalAlgorithm):
    @property
    def primary_binary(self):
        return "split-pipe"

    @property
    def name(self):
        return "split-pipe"

    @property
    def multi_core(self):
        return True  # fastqc has a threads option - and does not make use of it

    def get_version(self):
        env = os.environ.copy()
        if "LD_LIBRARY_PATH" in env:  # rpy2 likes to sneak this in, breaking e.g. STAR
            del env["LD_LIBRARY_PATH"]

        return (
            subprocess.check_output([self.primary_binary, "--version"], env=env)
            .decode("utf-8")
            .strip()
            .split()[-1]
        )

    def build_cmd(self, output_directory, ncores, arguments):
        res = [
            self.primary_binary,
        ] + arguments
        return res

    def prepare_index(self, genome):
        if not genome.name in _genome_job_cache:
            genome_name = genome.name + "_v1.0"
            output_dir = Path("cache/splitparse_genomes") / genome_name
            _genome_job_cache[genome.name] = self.run(
                output_dir,
                lambda: [
                    # fmt: off
                '-m', 'mkRef',
                '--genome_name', genome_name, # because if it's not in there, splitpipe fails
                #'--genome_dir', str(output_dir.absolute()),
                '--output_dir', str(output_dir.absolute()),
                '--fasta', str(genome.find_file('genome.fasta').absolute()),
                '--genes', str(genome.find_file('genes.gtf').absolute()),
                    # fmt: on
                ],
                additional_files_created=[
                    # output_dir / f"split-pipe_{self.get_version()}.log", # that has a ton of machine-dependend output.
                    output_dir / "all_genes.csv",
                    output_dir / "chrLength.txt",
                    output_dir / "chrNameLength.txt",
                    output_dir / "chrName.txt",
                    output_dir / "chrStart.txt",
                    output_dir / "exonGeTrInfo.tab",
                    output_dir / "exonInfo.tab",
                    output_dir / "exons.gtf.gz",
                    output_dir / "gene_info.json",
                    output_dir / "geneInfo.tab",
                    output_dir / "genes.gtf.gz",
                    output_dir / "Genome",
                    output_dir / "genome.fas.gz",
                    output_dir / "genomeParameters.txt",
                    output_dir / "SA",
                    output_dir / "SAindex",
                    output_dir / "sjdbInfo.txt",
                    output_dir / "sjdbList.fromGTF.out.tab",
                    output_dir / "sjdbList.out.tab",
                    output_dir / "transcriptInfo.tab",
                ],
            ).depends_on(genome.download())

        return _genome_job_cache[genome.name]

    def run_on_sublibrary(
        self, raw_sublibrary_sample, sample_definition, chemistry, genome
    ):
        # todo: add sample_definition sanity checker
        out_folder = Path("results/splitpipe") / raw_sublibrary_sample.name
        split_pipe_version = "1.1.2"

        additional_files_created = [
            # out_folder / f"split-pipe_{self.get_version()}.log", # that has a ton of machine dependend output
        ]
        input = []

        genome_job = self.prepare_index(genome)
        ncores = ppg.global_pipegraph.cores
        arguments = [
            # fmt: off
              '-m', 'all',
                '--chemistry', chemistry,
                '--nthreads', str(ncores),
                '--output_dir', str(out_folder.absolute()),
                '--genome_dir', str(genome_job.files[0].parent.absolute()),
            # fmt: on
        ]
        fastq_jobs = raw_sublibrary_sample.prepare_input_gzip()
        arguments += ["--fq1", str(fastq_jobs['R1'].files[0])]
        arguments += ["--fq2", str(fastq_jobs['R2'].files[0])]

        for sample, definition in sample_definition.items():
            additional_files_created.append(
                out_folder / f"{sample}_analysis_summary.html"
            )
            arguments += ["--sample", sample, definition]

        job = self.run(
            out_folder,
            arguments=arguments,
            additional_files_created=additional_files_created,
        )
        for file in input:
            job.depends_on_file(file)
        job.depends_on(genome_job)
        job.depends_on(fastq_jobs.values())
        return job
