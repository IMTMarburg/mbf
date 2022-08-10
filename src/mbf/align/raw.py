import pypipegraph as ppg
from pathlib import Path
from .strategies import build_fastq_strategy
from . import fastq2
from .exceptions import PairingError


class Sample:
    def __init__(
        self,
        sample_name,
        input_strategy,
        reverse_reads,
        fastq_processor=fastq2.Straight(),
        pairing="single",
        vid=None,
    ):
        """A sequenced sample, represented by one or more fastq files

        Paramaters
        ----------
            sample_name: str
                name of sample - must be unique
            input_strategy:  varied
                see build_fastq_strategy
            reverse_reads: bool
                whether to reverse the reads before processing
            fastq_processor: fastq2.*
                Preprocessing strategy
            pairing: 'auto', 'single', 'paired', 'only_first', 'only_second', 'paired_as_first'
                default: 'single'
                'auto' -> discover pairing from presence of R1/R2 files (-> 'single' or 'paired')
                'single' -> single end sequencing
                'paired' -> 'paired end' sequencing
                'only_first -> 'paired end' sequencing, but take only R1 reads
                'only_second' -> 'paired end' sequencing, but take only R2 reads
                'paired_as_single' -> treat each fragment as an independent read
                'paired_swap' -> 'paired end' sequencing, but swap R1/R2 reads
            vid: str
                sample identification number
        """
        self.name = sample_name
        ppg.assert_uniqueness_of_object(self)

        self.input_strategy = build_fastq_strategy(input_strategy)
        self.reverse_reads = reverse_reads
        self.fastq_processor = fastq_processor
        self.vid = vid
        accepted_pairing_values = (
            "auto",
            "single",
            "paired",
            "only_first",
            "only_second",
            "paired_as_single",
            "paired_swap",
        )
        if not pairing in accepted_pairing_values:
            raise ValueError(
                f"pairing was not in accepted values: {accepted_pairing_values}"
            )
        if pairing == "auto":
            if self.input_strategy.is_paired:
                pairing = "paired"
            else:
                pairing = "single"
        if pairing == "paired" and reverse_reads:
            raise ValueError(
                "pairing and inversing the reads does not play nice together, refusing."
            )
        self.pairing = pairing
        self.is_paired = self.pairing in ("paired", "paired_swap")
        self.cache_dir = (
            Path(ppg.util.global_pipegraph.cache_folder) / "lanes" / self.name
        )
        self.result_dir = Path("results") / "lanes" / self.name
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.result_dir.mkdir(parents=True, exist_ok=True)
        self.register_qc()

    def get_aligner_input_filenames(self):
        if self.is_paired:
            return (
                self.cache_dir / "input_R1_.fastq",
                self.cache_dir / "input_R2_.fastq",
            )
        else:
            return (self.cache_dir / "input.fastq",)

    def prepare_input(self):
        # input_strategy returns a list of
        # paired fastq files
        # ie. [('A_R1_.fastq1', 'A_R2.fastq', ...), ...]

        input_pairs = self.input_strategy()
        any_r2 = any([len(x) > 1 for x in input_pairs])
        # Single end - works from flat list
        if self.pairing == "single":
            if any_r2:
                raise PairingError(
                    f"{self.name}: paired end lane defined as single end - you need to change the pairing parameter"
                )
            input_filenames = [str(f[0]) for f in input_pairs]
        elif self.pairing == "paired_as_single":
            input_filenames = [str(f) for fl in input_pairs for f in fl]
        elif self.pairing == "only_first":
            input_filenames = [str(f[0]) for f in input_pairs]
        elif self.pairing == "only_second":
            input_filenames = [str(f[1]) for f in input_pairs]
        elif self.pairing == "paired":
            if not any_r2:
                raise PairingError(
                    f"Paired end lane, but no R2 reads found. Found files: {input_pairs}"
                )
            input_filenames = [
                (str(f[0]), str(f[1])) for f in input_pairs
            ] 
        elif self.pairing == "paired_swap":
            if not any_r2:
                raise PairingError(
                    f"Paired end lane, but no R2 reads found. Found files: {input_pairs}"
                )
            input_filenames = [
                (str(f[1]), str(f[0])) for f in input_pairs
            ]

        else:
            raise PairingError("unknown pairing")  # pragma: no cover
        if self.pairing in ["paired", "paired_swap"]:
            flat_input_filenames = [f for fl in input_pairs for f in fl]
        else:
            flat_input_filenames = input_filenames

        if hasattr(self.input_strategy, "dependencies"):
            deps = self.input_strategy.dependencies
        else:
            deps = [ppg.FileChecksumInvariant(f) for f in flat_input_filenames]
        output_filenames = self.get_aligner_input_filenames()

        if self.is_paired:
            if hasattr(self.fastq_processor, "generate_aligner_input_paired"):

                def prep_aligner_input():
                    import shutil

                    self.fastq_processor.generate_aligner_input_paired(
                        str(output_filenames[0]) + ".temp",
                        str(output_filenames[1]) + ".temp",
                        input_filenames,
                        self.reverse_reads,
                    )
                    shutil.move(str(output_filenames[0]) + ".temp", output_filenames[0])
                    shutil.move(str(output_filenames[1]) + ".temp", output_filenames[1])

                job = ppg.MultiTempFileGeneratingJob(
                    output_filenames, prep_aligner_input
                )
                job.depends_on(
                    self.fastq_processor.get_dependencies(
                        [str(x) for x in output_filenames]
                    )
                )
            else:

                def prep_aligner_input_r1():
                    import shutil

                    self.fastq_processor.generate_aligner_input(
                        str(output_filenames[0]) + ".temp",
                        [x[0] for x in input_filenames],
                        self.reverse_reads,
                    )
                    shutil.move(str(output_filenames[0]) + ".temp", output_filenames[0])

                def prep_aligner_input_r2():
                    import shutil

                    self.fastq_processor.generate_aligner_input(
                        str(output_filenames[1]) + ".temp",
                        [x[1] for x in input_filenames],
                        self.reverse_reads,
                    )
                    shutil.move(str(output_filenames[1]) + ".temp", output_filenames[1])

                jobR1 = ppg.TempFileGeneratingJob(
                    output_filenames[0], prep_aligner_input_r1
                )
                jobR2 = ppg.TempFileGeneratingJob(
                    output_filenames[1], prep_aligner_input_r2
                )

                jobR1.depends_on(
                    self.fastq_processor.get_dependencies(str(output_filenames[0]))
                )
                jobR2.depends_on(
                    self.fastq_processor.get_dependencies(str(output_filenames[1]))
                )
                job = ppg.JobList([jobR1, jobR2])
                # needed by downstream code.
                job.filenames = [output_filenames[0], output_filenames[1]]
        else:

            def prep_aligner_input(output_filename):
                import shutil

                self.fastq_processor.generate_aligner_input(
                    str(output_filename) + ".temp", input_filenames, self.reverse_reads
                )
                shutil.move(str(output_filename) + ".temp", output_filename)

            job = ppg.TempFileGeneratingJob(output_filenames[0], prep_aligner_input)
            job.depends_on(
                self.fastq_processor.get_dependencies(str(output_filenames[0]))
            )

        job.depends_on(
            deps,
            ppg.ParameterInvariant(
                self.name + "input_files",
                tuple(sorted(input_filenames))
                + (self.reverse_reads, self.fastq_processor.__class__.__name__),
            ),
        )
        return job

    def save_input(self):
        """Store the filtered input also in filename for later reference"""
        import gzip

        temp_job = self.prepare_input()
        output_dir = self.result_dir / "aligner_input"
        output_dir.mkdir(exist_ok=True)
        output_names = [output_dir / (Path(x).name + ".gz") for x in temp_job.filenames]
        pairs = zip(temp_job.filenames, output_names)

        def do_store():
            block_size = 10 * 1024 * 1024
            for input_filename, output_filename in pairs:
                op = open(input_filename, "rb")
                op_out = gzip.GzipFile(output_filename, "wb")
                f = op.read(block_size)
                while f:
                    op_out.write(f)
                    f = op.read(block_size)
                op_out.close()
                op.close()

        return ppg.MultiFileGeneratingJob(output_names, do_store).depends_on(temp_job)

    def align(self, aligner, genome, aligner_parameters, name=None):
        from .lanes import AlignedSample

        output_dir = (
            Path("results")
            / "aligned"
            / ("%s_%s" % (aligner.name, aligner.version))
            / genome.name
            / self.name
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        output_filename = output_dir / (self.name + ".bam")
        input_job = self.prepare_input()
        index_job = genome.build_index(aligner)

        alignment_job = aligner.align_job(
            input_job.filenames[0],
            input_job.filenames[1] if self.is_paired else None,
            index_job,
            output_filename,
            aligner_parameters if aligner_parameters else {},
        )
        alignment_job.depends_on(
            input_job,
            index_job,
            # ppg.ParameterInvariant(output_filename, aligner_parameters), # that's the aligner's job.
        )
        for j in alignment_job.prerequisites:
            if isinstance(j, ppg.ParameterInvariant):
                break
        else:
            raise ppg.JobContractError(
                "aligner (%s).align_job should have added a parameter invariant for aligner parameters"
                % aligner
            )
        return AlignedSample(
            f"{self.name if name is None else name}_{aligner.name}",
            alignment_job,
            genome,
            self.is_paired,
            self.vid,
            output_dir,
            aligner=aligner,
        )

    def register_qc(self):
        from mbf_qualitycontrol import qc_disabled

        if not qc_disabled():
            self.register_qc_fastqc()

    def register_qc_fastqc(self):
        from mbf_externals import FASTQC
        from mbf_qualitycontrol import register_qc

        a = FASTQC()
        output_dir = self.result_dir / "FASTQC"
        temp_job = self.prepare_input()
        if hasattr(temp_job, "filenames"):
            filenames = temp_job.filenames
        else:
            filenames = []
            for j in temp_job:  # is actually joblist
                filenames.extend(j.filenames)

        job = a.run(output_dir, filenames)
        return register_qc(job.depends_on(temp_job))
