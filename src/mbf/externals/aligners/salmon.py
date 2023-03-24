from ..externals import ExternalAlgorithm
import pypipegraph as ppg
from pathlib import Path
import hashlib


class Salmon(ExternalAlgorithm):
    def __init__(self):
        """@accepted_biotypes may be a set, or None to use all.

        While this is sorted into Aligners,
        it's not an aligner.
        use run_alevin_on_raw_lane(...)(

        """
        super().__init__()
        if accepted_biotypes is not None and not isinstance(accepted_biotypes, set):
            raise ValueError("@accepted_biotypes may be a set, or None to use all")
        self.accepted_biotypes = accepted_biotypes

    @property
    def name(self):
        return "Salmon"

    def primary_binary(self):
        return "salmon"

    @property
    def multi_core(self):
        return True

    def build_cmd(self, output_directory, ncores, arguments):
        return [self.primary_binary] + arguments + ["-p", str(ncores)]

    # def _aligner_build_cmd(self, output_dir, ncores, arguments):
    # return arguments + ["--runThreadN", str(ncores)]

    def get_index_version_range(self):  # pragma: no cover
        return None, None

    def get_genome_deps(self, genome):
        return [genome.job_transcripts()]

    def get_build_key(self):
        if self.accepted_biotypes is not None:
            return hashlib.md5(
                str(sorted(self.accepted_biotypes)).encode("utf-8")
            ).hexdigest()
        return hashlib.md5("None".encode("utf-8")).hexdigest()

    def build_index_from_genome(self, genome, output_fileprefix):
        import pysam

        output_dir = Path(output_fileprefix)
        output_key = self.get_build_key()
        output_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = Path("cache") / "salmon" / genome.name / output_key
        temp_dir.mkdir(parents=True, exist_ok=True)

        try:
            tf_cdna = open(temp_dir / "transcripts.fasta", "w")
            of_mapping = open(output_dir / "gene_transcript.mapping", "w")

            to_output = set()
            genes_out = set()

            for transcript in genome.transcripts.values():
                if (
                    self.accepted_biotypes is None
                    or transcript.biotype in self.accepted_biotypes
                ):
                    to_output.add(transcript.transcript_stable_id)
                    genes_out.add(transcript.gene_stable_id)
                    of_mapping.write(
                        f"{transcript.transcript_stable_id}\t{transcript.gene_stable_id}\n"
                    )

            seen = set()
            # iterating the fasta once is much faster than random accessing for each and every transcript
            with pysam.FastxFile(genome.find_file("cdna.fasta")) as f:
                for entry in f:
                    name = entry.name
                    if not name in to_output and "." in name:
                        name = name[: name.rfind(".")]
                    if name in to_output:
                        tf_cdna.write(f">{name}\n{entry.sequence}\n")
                        seen.add(name)
            if not seen:
                raise ValueError("non seen", entry.name)
            for transcript_stable_id in to_output.difference(seen):
                mrna = genome.transcripts[transcript_stable_id].mrna
                tf_cdna.write(f">{transcript_stable_id}\n{mrna}\n")

            with open(output_fileprefix / "mt.genes", "w") as op:
                for gene_stable_id in genome.df_genes.index[
                    genome.df_genes["name"].str.startswith("MT-")
                ]:
                    if gene_stable_id in genes_out:
                        op.write(gene_stable_id + "\n")
            with open(output_fileprefix / "rrna.genes", "w") as op:
                for gene_stable_id in genome.df_genes.index[
                    genome.df_genes["biotype"] == "rRNA"
                ]:
                    if gene_stable_id in genes_out:
                        op.write(gene_stable_id + "\n")

            tf_cdna.flush()
            of_mapping.close()
            run_func = self.get_run_func(
                output_fileprefix,
                [
                    "index",
                    "-i",
                    str((output_dir / "index").absolute()),
                    "-k",
                    "31",
                    "--transcripts",
                    str(Path(tf_cdna.name).absolute()),
                ],
            )
            run_func()

        finally:
            tf_cdna.close()
            # os.unlink(tf_cdna.name)
            of_mapping.close()

    def _run_alevin(
        self,
        output_path,
        fastqs_r1_r2,
        genome,
        method,
        dumpMtx=False,
        dumpFeatures=True,
    ):
        allowed_methods = (
            "chromium",
            "chromuimV3",
            "dropseq",
            "gemcode",
            "celseq",
            "celseq2",
            "quartzseq2",
        )
        if not method in allowed_methods:
            raise ValueError("method must be one of {allowed_methods}")
        output_path = Path(output_path)
        index_path = genome.build_index(self).output_path
        r1s = []
        r2s = []
        try:
            for r1, r2 in fastqs_r1_r2:
                r1s.append(str(r1))
                r2s.append(str(r2))
                if " " in str(r1) or " " in str(r2):
                    raise ValueError("Spaces in fastq names not allowed", r1, r2)
        except ValueError as e:
            if "not enough values" in str(e):
                raise ValueError(
                    "alevin currently works on paired end data like drop seq only?"
                )
            else:
                raise

        cmd = (
            [
                "alevin",
                "-lISR",  # inward, stranded, came from reverse strand... thats' right for dropseq
                "-1",
            ]
            + r1s
            + ["-2"]
            + r2s
            + [
                "--" + method,
                "-i",
                str(index_path / "index"),
                "-o",
                str(Path(output_path).absolute()),
                "--tgMap",
                str(index_path / "gene_transcript.mapping"),
                "--mrna",
                str(index_path / "mt.genes"),
                "--rrna",
                str(index_path / "rrna.genes"),
                "--maxNumBarcodes",
                "100000",  # just so you know what to patch
            ]
        )
        if dumpMtx:
            cmd.append("--dumpMtx")
        if dumpFeatures:
            cmd.append("--dumpFeatures")

        self.get_run_func(output_path, cmd)()

    def run_alevin_on_raw_lane(self, lane, genome, method):
        output = Path("results/alevin/") / lane.name

        def run_alevin():
            output.mkdir(exist_ok=True, parents=True)
            self._run_alevin(
                output, [lane.get_aligner_input_filenames()], genome, method
            )
            (output / "sentinel.txt").write_text("done")

        job = ppg.FileGeneratingJob(output / "sentinel.txt", run_alevin).depends_on(
            genome.build_index(self), lane.prepare_input()
        )
        job.depends_on_params(self.get_version())
        job.cores_needed = -1

        def run_qc():
            (output / "QC").mkdir(exist_ok=True)
            import rpy2.robjects as ro

            ro.r("library('alevinQC')")
            ro.r("alevinQCReport")(
                baseDir=str(output.absolute()),
                sampleId=lane.name,
                outputFile="alevinReport.html",
                outputFormat="html_document",
                outputDir=str((output / "QC").absolute()),
                forceOverwrite=True,
            )

        qc_job = ppg.FileGeneratingJob(
            output / "QC" / "alevinReport.html", run_qc
        ).depends_on(job)
        return job, qc_job

    def run_quant_on_raw_lane(
        self, lane, genome, libtype, options=None, gene_level=False
    ):
        output = Path(f"results/{self.name}/") / genome.name / lane.name

        def run_quant():
            output.mkdir(exist_ok=True, parents=True)
            self._run_quant(output, lane, genome, libtype, options, gene_level)
            (output / "sentinel.txt").write_text("done")

        job = ppg.FileGeneratingJob(output / "sentinel.txt", run_quant).depends_on(
            genome.build_index(self),
            lane.prepare_input(),
            ppg.FunctionInvariant("Salmon.run_quant", Salmon.run_quant),
        )
        job.depends_on_params(self.get_version())
        return job

    def _run_quant(
        self, outputpath, lane, genome, libtype, options=None, gene_level=False
    ):
        output_path = Path(outputpath)
        index_path = genome.build_index(self).output_path
        aligner_input = lane.get_aligner_input_filenames()
        cmd = ["quant", "-i", str(index_path / "index"), "-l", libtype]
        if gene_level:
            cmd.extend(["-g", str(index_path / "gene_transcript.mapping")])
        if lane.pairing == "single":
            cmd.extend(["-r", str(aligner_input[0])])
        elif lane.pairing == "paired":
            cmd.extend(["-1", str(aligner_input[0]), "-2", str(aligner_input[1])])
        elif lane.pairing == "only_first":
            cmd.extend(["-r", str(aligner_input[0])])
        elif lane.pairing == "only_second":
            cmd.extend(["-r", str(aligner_input[1])])
        elif lane.pairing == "paired_as_single":
            cmd.extend(["-r", str(aligner_input[0]) + " " + str(aligner_input[1])])
        else:
            raise ValueError("Could not understand pairing mode: %s" % lane.pairing)
        cmd.extend(
            [
                "--validateMappings",
                "-o",
                str(output_path),
            ]
        )
        if options is not None:
            for key in options:
                if options[key] is not None:
                    cmd.extend([key, options[key]])
                else:
                    cmd.append(key)
        self.get_run_func(output_path, cmd)()

    def get_index_filenames(self):
        [
            "index/ctable.bin",
            "index/ctg_offsets.bin",
            "index/duplicate_clusters.tsv",
            "index/eqtable.bin",
            "index/info.json",
            "index/mphf.bin",
            "index/pos.bin",
            "index/pre_indexing.log",
            "index/rank.bin",
            "index/refAccumLengths.bin",
            "index/ref_indexing.log",
            "index/reflengths.bin",
            "index/refseq.bin",
            "index/seq.bin",
            "index/versionInfo.json",
        ]
