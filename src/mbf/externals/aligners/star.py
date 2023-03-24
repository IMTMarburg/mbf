from .base import Aligner
import pypipegraph as ppg
from pathlib import Path
import subprocess


class STAR(Aligner):
    @property
    def name(self):
        return "STAR"

    @property
    def primary_binary(self):
        return "STAR"

    @property
    def multi_core(self):
        return True

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        return arguments + ["--runThreadN", str(ncores)]

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_job,
        output_bam_filename,
        parameters,
    ):
        def build_cmd():
            """must be delayed for the index job..."""
            if hasattr(index_job, "target_folder"):  # ppg2 sharedmultifilegenjob
                index_path = index_job.target_folder
            elif hasattr(index_job, "output_path"):  # ppg1 PrebuildJob
                index_path = index_job.output_path
            else:
                index_path = Path(index_job.files[0]).parent

            cmd = [
                "FROM_ALIGNER",
                "STAR",
                "--genomeDir",
                str(Path(index_path).absolute()),
                "--genomeLoad",
                "NoSharedMemory",
                "--readFilesIn",
            ]
            if "," in str(input_fastq) or (
                paired_end_filename and "," in str(paired_end_filename)
            ):  # pragma: no cover
                raise ValueError("STAR does not handle fastq filenames with a comma")
            if paired_end_filename:
                cmd.extend(
                    [
                        '"%s"' % Path(paired_end_filename).absolute(),
                        '"%s"' % Path(input_fastq).absolute(),
                    ]
                )
            else:
                cmd.extend([Path(input_fastq).absolute()])
            cmd.extend(["--outSAMtype", "BAM", "SortedByCoordinate"])
            for k, v in parameters.items():
                cmd.append(k)
                cmd.append(str(v))
            return cmd

        additional_files_created = [output_bam_filename]
        if parameters.get("--outReadsUnmapped", False) == "Fastx":
            additional_files_created.append(
                Path(output_bam_filename).parent / "Unmapped.out.mate1"
            )
            if paired_end_filename:
                additional_files_created.append(
                    Path(output_bam_filename).parent / "Unmapped.out.mate2"
                )

        def rename_after_alignment():
            ob = Path(output_bam_filename)
            (ob.parent / "Aligned.sortedByCoord.out.bam").rename(ob.parent / ob.name)

        job = self.run(
            Path(output_bam_filename).parent,
            build_cmd,
            cwd=Path(output_bam_filename).parent,
            call_afterwards=rename_after_alignment,
            additional_files_created=additional_files_created,
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items())),
            index_job,
        )
        return job

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        if isinstance(fasta_files, (str, Path)):
            fasta_files = [fasta_files]
        if len(fasta_files) > 1:
            raise ValueError("STAR can only build from a single fasta")
        # if gtf_input_filename is None:
        #     raise ValueError(
        #         "STAR needs a gtf input file to calculate splice junctions"
        #     )
        cmd = [
            "FROM_ALIGNER",
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            Path(output_fileprefix).absolute(),
            "--genomeFastaFiles",
            Path(fasta_files[0]).absolute(),
        ]
        if (
            gtf_input_filename is not None
            and Path(gtf_input_filename).stat().st_size > 0
        ):
            cmd.extend(
                [
                    "--sjdbGTFfile",
                    Path(gtf_input_filename).absolute(),
                    "--sjdbOverhang",
                    "100",
                ]
            )
        if Path(fasta_files[0]).stat().st_size < 5e4:  # small genome,
            # adjust annoying parameter
            # no clue why star doesn't do this automatically
            import pysam
            import math

            # straight from the manual
            total = 0
            for entry in pysam.FastxFile(fasta_files[0]):
                total += len(entry.sequence)
            print("total genome length", total)
            genomeSAindexNbases = min(14, math.floor(math.log2(total) / 2) - 1)
            cmd.extend(["--genomeSAindexNbases", str(genomeSAindexNbases)])

        def fix_up_if_gtf_was_empty():
            # if there was no gtf, or the GTF was empty (e.g. phix)
            # we don't have these files.
            # so we fake them
            if (
                gtf_input_filename is None
                or Path(gtf_input_filename).stat().st_size == 0
            ):
                for fn in [
                    "exonGeTrInfo.tab",
                    "exonInfo.tab",
                    "geneInfo.tab",
                    "sjdbInfo.txt",
                    "sjdbList.fromGTF.out.tab",
                    "sjdbList.out.tab",
                    "transcriptInfo.tab",
                ]:
                    (output_fileprefix / fn).write_text("")

        return self.get_run_func(
            output_fileprefix,
            cmd,
            cwd=output_fileprefix,
            call_afterwards=fix_up_if_gtf_was_empty,
        )

    def get_version(self):
        return subprocess.check_output(["STAR", "--version"]).decode("utf-8").strip()

    def get_alignment_stats(self, output_bam_filename):
        target = Path(output_bam_filename).parent / "Log.final.out"
        if not target.exists():  # pragma: no cover
            return {"No data found": 1}
        else:
            lines = target.read_text().split("\n")
            lines = [x.split(" |", 1) for x in lines if " |" in x]
            lookup = {x[0].strip(): x[1].strip() for x in lines}
            result = {}
            for k in [
                "Number of reads mapped to too many loci",
                "Uniquely mapped reads number",
                "Number of reads mapped to multiple loci",
            ]:
                result[k] = int(lookup[k])
            result["Unmapped"] = int(lookup["Number of input reads"]) - sum(
                result.values()
            )
            return result

    def get_index_filenames(self):
        res = [
            "chrLength.txt",
            "chrNameLength.txt",
            "chrName.txt",
            "chrStart.txt",
            "exonGeTrInfo.tab",
            "exonInfo.tab",
            "geneInfo.tab",
            "Genome",
            "genomeParameters.txt",
            "SA",
            "SAindex",
            "sjdbInfo.txt",
            "sjdbList.fromGTF.out.tab",
            "sjdbList.out.tab",
            "transcriptInfo.tab",
        ]
        return res


class STARSolo(STAR):
    @property
    def name(self):
        return "STARSolo"

    def _parse_cellwhitelist_job(self, job):
        if isinstance(job, list):
            return [self._parse_cellwhitelist_job(x) for x in job]
        elif isinstance(job, (str, Path)):
            return ppg.FileInvariant(job)
        if not isinstance(
            job,
            (ppg.FileInvariant, ppg.FileGeneratingJob, ppg.MultiFileGeneratingJob),
        ):
            if hasattr(ppg, "is_ppg2"):
                import pypipegraph2 as ppg2

                if isinstance(
                    job,
                    (
                        ppg2.FileInvariant,
                        ppg2.FileGeneratingJob,
                        ppg2.MultiFileGeneratingJob,
                    ),
                ):
                    return job
            raise ValueError(
                "job must be a ppg.FileGeneratingJob or FileChecksumInvariant"
                "was %s" % (type(job))
            )
        else:
            return job

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_job,
        output_bam_filename,
        parameters,
    ):
        if "cell_barcode_whitelist" in parameters:
            cell_whitelist_job = self._parse_cellwhitelist_job(
                parameters["cell_barcode_whitelist"]
            )
            if not isinstance(cell_whitelist_job, list):
                cell_whitelist_job = [cell_whitelist_job]
            del parameters["cell_barcode_whitelist"]
        else:
            cell_whitelist_job = None

        soloType = parameters["soloType"]
        del parameters["soloType"]
        allowed_solotypes = ("CB_UMI_Simple", "CB_UMI_Complex")
        if not soloType in allowed_solotypes:
            raise ValueError(
                "unsupported solo type", soloType, "allowed", allowed_solotypes
            )

        gtf_job = None
        if "sjdbGTFfile" in parameters:
            if isinstance(parameters['sjdbGTFfile'], ppg.Job):
                gtf_job = parameters['sjdbGTFfile']
                parameters['sjdbGTFfile'] = Path(gtf_job.job_id).absolute()
            elif isinstance(parameters['sjdbGTFfile'], (Path, str)):
                gtf_job = ppg.FileInvariant(parameters['sjdbGTFfile'])
                parameters['sjdbGTFfile'] = Path(gtf_job.job_id).absolute()

        def build_cmd():
            """must be delayed for the index job..."""
            if hasattr(index_job, "target_folder"):  # ppg2 sharedmultifilegenjob
                index_path = index_job.target_folder
            elif hasattr(index_job, "output_path"):  # ppg1 PrebuildJob
                index_path = index_job.output_path
            else:
                index_path = Path(index_job.files[0]).parent

            cmd = [
                "FROM_ALIGNER",
                "STAR",
                "--genomeDir",
                str(Path(index_path).absolute()),
                "--genomeLoad",
                "NoSharedMemory",
                "--soloType",
                soloType,
            ]

            if cell_whitelist_job is not None:
                cmd.append(
                    "--soloCBwhitelist",
                )
                for j in cell_whitelist_job:
                    cmd.append(str(j.files[0].absolute()))
            else:
                "--soloCBwhitelist=None",

            cmd.append(
                "--readFilesIn",
            )
            if "," in str(input_fastq) or (
                paired_end_filename and "," in str(paired_end_filename)
            ):  # pragma: no cover
                raise ValueError("STAR does not handle fastq filenames with a comma")
            if paired_end_filename:
                cmd.extend(
                    [
                        # not sure. the doc says read2 then read1,
                        # cdna, then barcode
                        # but our read1 is the barcode+umi, so
                        '"%s"' % Path(paired_end_filename).absolute(),
                        '"%s"' % Path(input_fastq).absolute(),
                    ]
                )
            else:
                raise ValueError("expected paired end data")

            cmd.extend(
                [
                    "--outSAMattributes",
                    "NH",
                    "HI",
                    "nM",
                    "AS",
                    "CR",  # raw (uncorrected) cellbarcode/umi
                    "UR",
                    "CB",  # corrected barcode and umi
                    "UB",
                    "GX",
                    "GN",  # gene id
                    "sS",  # combined cellbarcode and umi
                    "sQ",  #
                    "sM",
                    "--outSAMtype",
                    "BAM",
                    "SortedByCoordinate",
                ]
            )

            for k, v in parameters.items():
                cmd.append(k if k.startswith("-") else "--" + k)
                if isinstance(v, list):
                    for vx in v:
                        cmd.append(str(vx))
                else:
                    cmd.append(str(v))
            return cmd

        output_directory = Path(output_bam_filename).parent

        additional_files_created = [
            output_bam_filename,
            output_directory / "Solo.out/Barcodes.stats",
            output_directory / "Solo.out/Gene/Features.stats",
            output_directory / "Solo.out/Gene/Summary.csv",
            output_directory / "Solo.out/Gene/UMIperCellSorted.txt",
            output_directory / "Solo.out/Gene/filtered/barcodes.tsv.gz",
            output_directory / "Solo.out/Gene/filtered/features.tsv.gz",
            output_directory / "Solo.out/Gene/filtered/matrix.mtx.gz",
            output_directory / "Solo.out/Gene/raw/barcodes.tsv.gz",
            output_directory / "Solo.out/Gene/raw/features.tsv.gz",
            output_directory / "Solo.out/Gene/raw/matrix.mtx.gz",
        ]  # todo: double check if that's all.

        if parameters.get("--outReadsUnmapped", False) == "Fastx":
            additional_files_created.append(
                Path(output_bam_filename).parent / "Unmapped.out.mate1"
            )
            if paired_end_filename:
                additional_files_created.append(
                    Path(output_bam_filename).parent / "Unmapped.out.mate2"
                )

        def rename_after_alignment():
            ob = Path(output_bam_filename)
            (ob.parent / "Aligned.sortedByCoord.out.bam").rename(ob.parent / ob.name)
            subprocess.check_call(
                [
                    "gzip",
                    output_directory / "Solo.out/Gene/raw/matrix.mtx",
                ]
            )
            subprocess.check_call(
                [
                    "gzip",
                    output_directory / "Solo.out/Gene/raw/features.tsv",
                ]
            )
            subprocess.check_call(
                [
                    "gzip",
                    output_directory / "Solo.out/Gene/raw/barcodes.tsv",
                ]
            )
            subprocess.check_call(
                [
                    "gzip",
                    output_directory / "Solo.out/Gene/filtered/matrix.mtx",
                ]
            )
            subprocess.check_call(
                [
                    "gzip",
                    output_directory / "Solo.out/Gene/filtered/features.tsv",
                ]
            )
            subprocess.check_call(
                [
                    "gzip",
                    output_directory / "Solo.out/Gene/filtered/barcodes.tsv",
                ]
            )

        job = self.run(
            Path(output_bam_filename).parent,
            build_cmd,
            cwd=Path(output_bam_filename).parent,
            call_afterwards=rename_after_alignment,
            additional_files_created=additional_files_created,
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items())),
            index_job,
        )
        if cell_whitelist_job is not None:
            job.depends_on(cell_whitelist_job)
        if gtf_job is not None:
            job.depends_on(gtf_job)
        return job
