# forwards for compatibility with old chipseq code

from .raw import Sample
import pypipegraph as ppg
import pysam
from pathlib import Path
import pandas as pd
import collections
from dppd import dppd
import dppd_plotnine  # noqa:F401 -
from mbf_qualitycontrol import register_qc, QCCollectingJob, qc_disabled

dp, X = dppd()


class _ChromosomeMangledSamFile(pysam.Samfile):
    """Wraps a samfile so that it understands targets that don't quite
    have the right name (eg. missing 'chr', additional 'chr' in front
    of chromosomes, etc.
    Usage:
    b = _ChromosomeMangledSamFile(("my.bam",'rb')
    b.chromosome_mangler = lambda x: 'chr' + x

    the chromosome mangler may return False, which means 'that's a legitimate region, but missing here - just act as if there are no reads on it'
    """

    def parse_region(
        self,
        contig=None,
        start=None,
        stop=None,
        region=None,
        tid=None,
        reference=None,
        end=None,
    ):
        if not hasattr(self, "chromosome_mangler"):
            raise ValueError(
                "You need to set a .chromosome_mangler on ChromosomeMangledSamFiles. Sorry about that - could not extend c__init__"
            )
        if reference:  # old name...
            contig = reference
            reference = None
        if end:
            stop = end
            end = None
        if contig:
            contig = self.chromosome_mangler(contig)
            print("new reference", contig)
            if contig is False:
                raise ValueError(
                    "Chromosome mangler for %s returned a region not in the file"
                    % (self,)
                )
        return pysam.Samfile.parse_region(
            self,
            contig=contig,
            start=start,
            stop=stop,
            region=region,
            tid=tid,
            reference=None,
            end=None,
        )


class _BamDerived:
    def _parse_alignment_job_input(self, alignment_job):
        if isinstance(alignment_job, (str, Path)):
            alignment_job = ppg.FileInvariant(alignment_job)
        if not isinstance(
            alignment_job,
            (ppg.FileInvariant, ppg.FileGeneratingJob, ppg.MultiFileGeneratingJob),
        ):
            raise ValueError(
                "alignment_job must be a ppg.FileGeneratingJob or FileChecksumInvariant"
                "was %s" % (type(alignment_job))
            )
        bam_name = None
        bai_name = None
        for fn in alignment_job.filenames:
            if str(fn).endswith(".bam"):
                if bam_name is None:
                    bam_name = str(fn)
                else:
                    raise ValueError(
                        "Job passed to AlignedSample had multiple .bam filenames"
                    )
            elif str(fn).endswith(".bai"):
                if bai_name is None:
                    index_fn = str(fn)
                    bai_name = index_fn
                else:
                    raise ValueError(
                        "Job passed to AlignedSample had multiple .bai filenames"
                    )

        if bam_name is None:
            raise ValueError("Job passed to AlignedSample had no .bam filenames")

        if isinstance(alignment_job, ppg.MultiFileGeneratingJob):
            if bai_name is None:
                index_fn = bam_name + ".bai"
                index_job = ppg.FileGeneratingJob(
                    index_fn, self._index(bam_name, index_fn)
                )
                index_job.depends_on(alignment_job)

            else:
                index_fn = bai_name
                index_job = alignment_job

        elif isinstance(alignment_job, ppg.FileGeneratingJob):
            index_fn = bam_name + ".bai"
            index_job = ppg.FileGeneratingJob(index_fn, self._index(bam_name, index_fn))
            index_job.depends_on(alignment_job)
        elif isinstance(alignment_job, ppg.FileInvariant):
            index_fn = bam_name + ".bai"
            if Path(index_fn).exists():
                index_job = ppg.FileInvariant(index_fn)
            else:
                cache_dir = Path(ppg.util.global_pipegraph.cache_folder) / "bam_indices"
                cache_dir.mkdir(exist_ok=True, parents=True)
                index_fn = cache_dir / (self.name + "_" + Path(bam_name).name + ".bai")
                index_job = ppg.FileGeneratingJob(
                    index_fn, self._index(bam_name, index_fn)
                )
                index_job.depends_on(alignment_job)
        else:
            raise NotImplementedError("Should not happe / covered by earlier if")
        return alignment_job, index_job, Path(bam_name), Path(index_fn)

    def _index(self, input_fn, output_fn):
        def do_index():
            pysam.index(str(Path(input_fn).absolute()), str(Path(output_fn).absolute()))

        return do_index

    def __hash__(self):
        return hash(self.__class__.__name__ + self.name)

    def load(self):
        return self.alignment_job, self.index_job, self.genome.download_genome

    def get_bam(self):
        import multiprocessing

        mapper = getattr(self, "chromosome_mapper", None)
        if mapper is not None:
            r = _ChromosomeMangledSamFile(
                self.bam_filename,
                index_filename=str(self.index_filename),
                threads=multiprocessing.cpu_count(),
            )
            r.chromosome_mangler = self.chromosome_mapper
        else:
            r = pysam.Samfile(
                self.bam_filename,
                index_filename=str(self.index_filename),
                threads=multiprocessing.cpu_count(),
            )
        return r

    def get_bam_names(self):
        """Retrieve the bam filename and index name as strings"""
        return (str(self.bam_filename), str(self.index_filename))


class AlignedSample(_BamDerived):
    def __init__(
        self,
        name,
        alignment_job,
        genome,
        is_paired,
        vid,
        result_dir=None,
        aligner=None,
        chromosome_mapper=None,
    ):
        """
        Create an aligned sample from a BAM producing job.
        See Sample.align()

        Parameters:
            alignment_job: FileGeneratingJob, FileInvariant, str, pathlib.Path
                Where does the BAM come from?
                str and Path get's converted into a FileInvariant
            genome:
                an mbf_genomes.* Genome
            is_paired: bool
                whether this is a paired end sequencing run
            vid: str
                a unique, external sample management id
            chromosome_mapper: Option[function]
                A function mapping pipeline chromosomes to those used in the BAM
        """

        self.name = name
        ppg.util.assert_uniqueness_of_object(self)
        (
            self.alignment_job,
            self.index_job,
            bam_name,
            index_fn,
        ) = self._parse_alignment_job_input(alignment_job)
        self.result_dir = (
            Path(result_dir)
            if result_dir
            else (Path("results") / "aligned" / self.name)
        )
        self.result_dir.mkdir(exist_ok=True, parents=True)
        self.genome = genome
        self.is_paired = is_paired
        self.vid = vid
        self.bam_filename = bam_name
        self.index_filename = index_fn
        self.aligner = aligner
        self.register_qc()
        self.chromosome_mapper = chromosome_mapper

    def get_unique_aligned_bam(self):
        """Deprecated compability with older pipeline"""
        return self.get_bam()

    def _parse_idxstat(self):
        by_chr = self.get_bam().get_index_statistics()
        mapped = 0
        unmapped = 0
        for record in by_chr:
            mapped += record.mapped
            unmapped += record.unmapped
        return mapped, unmapped

    def mapped_reads(self):
        """How many mapped entrys are in the bam?"""
        return self._parse_idxstat()[0]

    def unmapped_reads(self):
        """How many unmapped entrys are in the bam?"""
        return self._parse_idxstat()[1]

    def post_process(self, post_processor, new_name=None, result_dir=None):
        """Postprocess this lane using a  mbf_align.postprocess.*
        Ie. Turn a lane into a 'converted' lane.

        """
        if new_name is None:
            new_name = self.name + "_" + post_processor.name
        if result_dir is None:
            result_dir = (
                self.result_dir
                / ".."
                / post_processor.result_folder_name
                / self.result_dir.name
            )
        result_dir = Path(result_dir)
        result_dir.mkdir(exist_ok=True, parents=True)
        bam_filename = result_dir / (new_name + ".bam")

        def inner(output_filename):
            post_processor.process(
                Path(self.get_bam_names()[0]), Path(output_filename), result_dir
            )

        alignment_job = ppg.FileGeneratingJob(bam_filename, inner)
        alignment_job.depends_on(
            self.load(),
            post_processor.get_dependencies(),
            ppg.ParameterInvariant(alignment_job.job_id, post_processor.get_parameters()),
        )
        vid = post_processor.get_vid(self.vid)

        new_lane = AlignedSample(
            new_name,
            alignment_job,
            self.genome,
            self.is_paired,
            vid,
            result_dir=result_dir,
        )

        new_lane.post_processor_jobs = post_processor.further_jobs(new_lane, self)
        new_lane.parent = self
        new_lane.post_processor_qc_jobs = post_processor.register_qc(new_lane)
        return new_lane

    def to_fastq(self, output_filename, as_temp_file=False):
        """Convert a (single end) bam back into a fastq"""
        if self.is_paired:
            raise ValueError(
                "No support for -> fastq for paired end files at the the moment"
            )

        def convert():
            import mbf_bam

            mbf_bam.bam_to_fastq(output_filename, self.get_bam_names()[0])

        if as_temp_file:
            cls = ppg.TempFileGeneratingJob
        else:
            cls = ppg.FileGeneratingJob
        return cls(output_filename, convert).depends_on(self.load())

    def get_alignment_stats(self):
        if self.aligner is not None and hasattr(self.aligner, "get_alignment_stats"):
            return self.aligner.get_alignment_stats(Path(self.bam_filename))
        else:
            with self.get_bam() as f:
                return {"Mapped": f.mapped, "Unmapped": f.unmapped}

    def register_qc(self):
        if not qc_disabled():
            self.register_qc_complexity()
            self.register_qc_gene_strandedness()
            self.register_qc_biotypes()
            self.register_qc_alignment_stats()
            self.register_qc_subchromosomal()
            self.register_qc_splicing()

    def register_qc_complexity(self):

        output_filename = self.result_dir / f"{self.name}_complexity.png"

        def calc():
            import mbf_bam

            counts = mbf_bam.calculate_duplicate_distribution(
                str(self.bam_filename), str(self.index_filename)
            )
            return pd.DataFrame(
                {
                    "source": self.name,
                    "Repetition count": list(counts.keys()),
                    "Count": list(counts.values()),
                }
            )

        def plot(df):
            import numpy as np

            unique_count = df["Count"].sum()
            total_count = (df["Count"] * df["Repetition count"]).sum()
            pcb = float(unique_count) / total_count
            if pcb >= 0.9:  # pragma: no cover
                severity = "none"
            elif pcb >= 0.8:  # pragma: no cover
                severity = "mild"
            elif pcb >= 0.5:  # pragma: no cover
                severity = "moderate"
            else:
                severity = "severe"
            title = (
                "Genomic positions with repetition count reads\nTotal read count: %i\nPCR Bottleneck coefficient: %.2f (%s)"
                % (total_count, pcb, severity)
            )
            return (
                dp(df)
                .p9()
                .theme_bw()
                .add_point("Repetition count", "Count")
                .add_line("Repetition count", "Count")
                .scale_y_continuous(
                    trans="log2",
                    breaks=[2 ** x for x in range(1, 24)],
                    labels=lambda x: ["2^%0.f" % np.log(xs) for xs in x],
                )
                .title(title)
                .pd
            )

        return register_qc(
            ppg.PlotJob(output_filename, calc, plot)
            .depends_on(self.load())
            .use_cores(-1)
        )

    def register_qc_gene_strandedness(self):  # noqa: C901
        from mbf_genomics.genes.anno_tag_counts import _IntervalStrategy

        class IntervalStrategyExonIntronClassification(_IntervalStrategy):
            """For QC purposes, defines all intron/exon intervals tagged
            with nothing but intron/exon

            See mbf_align.lanes.AlignedLane.register_qc_gene_strandedness

            """

            def _get_interval_tuples_by_chr(self, genome):
                from mbf_nested_intervals import IntervalSet

                coll = {chr: [] for chr in genome.get_chromosome_lengths()}
                for g in genome.genes.values():
                    exons = g.exons_overlapping
                    if len(exons[0]) == 0:  # pragma: no cover
                        exons = g.exons_merged
                    for start, stop in zip(*exons):
                        coll[g.chr].append(
                            (start, stop, 0b0101 if g.strand == 1 else 0b0110)
                        )
                    for start, stop in zip(*g.introns_strict):
                        coll[g.chr].append(
                            (start, stop, 0b1001 if g.strand == 1 else 0b1010)
                        )
                result = {}
                for chr, tups in coll.items():
                    iset = IntervalSet.from_tuples_with_id(tups)
                    # iset = iset.merge_split()
                    iset = iset.merge_hull()
                    if iset.any_overlapping():
                        raise NotImplementedError("Should not be reached")
                    result[chr] = []
                    for start, stop, ids in iset.to_tuples_with_id():
                        ids = set(ids)
                        if len(ids) == 1:
                            id = list(ids)[0]
                            if id == 0b0101:
                                tag = "exon"
                                strand = +1
                            elif id == 0b0110:
                                tag = "exon"
                                strand = -1
                            elif id == 0b1001:
                                tag = "intron"
                                strand = +1
                            elif id == 0b1010:
                                tag = "intron"
                                strand = -1
                            else:  # pragma: no cover
                                raise NotImplementedError("Should not be reached")
                        else:
                            down = 0
                            for i in ids:
                                down |= i
                            if down & 0b1100 == 0b1100:
                                tag = "both"
                            elif down & 0b0100 == 0b0100:
                                tag = "exon"
                            else:  # pragma: no cover  haven't observed this case in the wild yet.
                                tag = "intron"  # pragma: no cover  # pragma: no cover  # pragma: no cover  haven't observed this case in the wild yet.
                            if down & 0b11 == 0b11:
                                tag += "_undecidable"
                                strand = (
                                    1  # doesn't matter, but must be one or the other
                                )
                            elif down & 0b01:
                                strand = 1
                            else:
                                strand -= 1

                        result[chr].append((tag, strand, [start], [stop]))
                return result

        output_filename = self.result_dir / f"{self.name}_strandedness.png"

        def calc():
            from mbf_genomics.genes.anno_tag_counts import IntervalStrategyGene
            from mbf_bam import count_reads_stranded

            interval_strategy = IntervalStrategyExonIntronClassification()
            intervals = interval_strategy._get_interval_tuples_by_chr(self.genome)

            bam_filename, bam_index_name = self.get_bam_names()
            forward, reverse = count_reads_stranded(
                bam_filename,
                bam_index_name,
                intervals,
                IntervalStrategyGene()._get_interval_tuples_by_chr(self.genome),
                each_read_counts_once=True,
            )
            result = {"what": [], "count": [], "sample": self.name}
            for k in forward.keys() | reverse.keys():
                if k.endswith("_undecidable"):
                    result["what"].append(k)
                    result["count"].append(forward.get(k, 0) + reverse.get(k, 0))
                elif not k.startswith("_"):
                    result["what"].append(k + "_correct")
                    result["count"].append(forward.get(k, 0))
                    result["what"].append(k + "_reversed")
                    result["count"].append(reverse.get(k, 0))
                elif k == "_outside":
                    result["what"].append("outside")
                    result["count"].append(forward.get(k, 0))

            return pd.DataFrame(result)

        def plot(df):
            return (
                dp(df)
                .mutate(
                    what=pd.Categorical(
                        df["what"],
                        [
                            "exon_correct",
                            "exon_reversed",
                            "exon_undecidable",
                            "intron_correct",
                            "intron_reversed",
                            "intron_undecidable",
                            "both_correct",
                            "both_reversed",
                            "both_undecidable",
                            "outside",
                        ],
                    )
                )
                .p9()
                .add_bar("sample", "count", fill="what", position="dodge")
                .scale_y_continuous(labels=lambda xs: ["%.2g" % x for x in xs])
                .turn_x_axis_labels()
                .pd
            )

        return register_qc(
            ppg.PlotJob(output_filename, calc, plot)
            .depends_on(self.load())
            .use_cores(-1)
        )

    def register_qc_biotypes(self):
        output_filename = self.result_dir / f"{self.name}_reads_per_biotype.png"

        from mbf_genomics.genes import Genes
        from mbf_genomics.genes.anno_tag_counts import GeneUnstranded

        genes = Genes(self.genome)
        anno = GeneUnstranded(self)

        def plot(output_filename):
            print(genes.df.columns)
            return (
                dp(genes.df)
                .groupby("biotype")
                .summarize((anno.columns[0], lambda x: x.sum(), "read count"))
                .mutate(sample=self.name)
                .p9()
                .theme_bw()
                .annotation_stripes()
                .add_bar("biotype", "read count", stat="identity")
                .scale_y_continuous(labels=lambda xs: ["%.2g" % x for x in xs])
                # .turn_x_axis_labels()
                .coord_flip()
                .title(self.name)
                .render(
                    output_filename,
                    width=6,
                    height=2 + len(genes.df.biotype.unique()) * 0.25,
                )
            )

        return register_qc(
            ppg.FileGeneratingJob(output_filename, plot).depends_on(
                genes.add_annotator(anno)
            )
        )

    def register_qc_alignment_stats(self):
        output_filename = self.result_dir / ".." / "alignment_statistics.png"

        def calc_and_plot(output_filename, lanes):
            parts = []
            for lane in lanes:
                p = lane.get_alignment_stats()
                parts.append(
                    pd.DataFrame(
                        {
                            "what": list(p.keys()),
                            "count": list(p.values()),
                            "sample": lane.name,
                        }
                    )
                )
            df = pd.concat(parts)
            order = sorted(df["what"].unique())
            umrn = "Uniquely mapped reads number"
            if umrn in order:
                order = [x for x in order if x != umrn] + [umrn]
            df.to_pickle(Path(output_filename).with_suffix('.pickle'))
            return (
                dp(df)
                .categorize("what", order)
                .p9()
                .theme_bw()
                .annotation_stripes()
                .add_bar(
                    "sample", "count", position="stack", stat="identity"
                )
                .facet_wrap("what", ncol=1)
                .title(lanes[0].genome.name)
                .turn_x_axis_labels()
                .scale_y_continuous(labels=lambda xs: ["%.2g" % x for x in xs])
                .render_args(width=len(parts) * 0.2 + 1, height=5, limitsize=False)
                .render(output_filename)
            )

        return register_qc(
            QCCollectingJob(output_filename, calc_and_plot)
            .depends_on(self.load())
            .add(self)
        )  # since everybody says self.load, we get them all

    def register_qc_subchromosomal(self):
        """Subchromosom distribution plot - good to detect amplified regions
        or ancient virus awakening"""
        import mbf_genomics

        output_filename = (
            self.result_dir / f"{self.name}_subchromosomal_distribution.png"
        )

        class IntervalStrategyWindows(
            mbf_genomics.genes.anno_tag_counts._IntervalStrategy
        ):
            """For QC purposes, spawn all chromosomes with
            windows of the definied size

            See mbf_align.lanes.AlignedLane.register_qc_subchromosomal

            """

            def __init__(self, window_size):
                self.window_size = window_size

            def _get_interval_tuples_by_chr(self, genome):
                result = {}
                for chr, length in genome.get_chromosome_lengths().items():
                    result[chr] = []
                    for ii in range(0, length, self.window_size):
                        result[chr].append(
                            ("%s_%i" % (chr, ii), 0, [ii], [ii + self.window_size])
                        )
                return result

        def calc():
            from mbf_bam import count_reads_unstranded

            interval_strategy = IntervalStrategyWindows(250_000)
            intervals = interval_strategy._get_interval_tuples_by_chr(self.genome)

            bam_filename, bam_index_name = self.get_bam_names()
            counts = count_reads_unstranded(
                bam_filename,
                bam_index_name,
                intervals,
                intervals,
                each_read_counts_once=True,
            )
            true_chromosomes = set(self.genome.get_true_chromosomes())
            result = {"chr": [], "window": [], "count": []}
            for key, count in counts.items():
                if not key.startswith("_"):
                    # must handle both 2R_1234
                    # and Unmapped_scaffold_29_D1705_1234
                    *c, window = key.split("_")
                    chr = "_".join(c)
                    if chr in true_chromosomes:  # pragma: no branch
                        window = int(window)
                        result["chr"].append(chr)
                        result["window"].append(window)
                        result["count"].append(count)
            return pd.DataFrame(result)

        def plot(df):
            import natsort

            df[
                "count"
            ] += 1  # so we don't crash in the log scale if all values are 0 for a chr
            return (
                dp(df)
                .categorize("chr", natsort.natsorted(X["chr"].unique()))
                .p9()
                .theme_bw()
                .add_line("window", "count", _alpha=0.3)
                .scale_y_log10()
                .facet_wrap("chr", scales="free", ncol=1)
                .hide_x_axis_labels()
                .title(self.name)
                .render_args(
                    width=6, height=2 + len(df["chr"].unique()) * 1, limitsize=False
                )
                .pd
            )

        return register_qc(
            ppg.PlotJob(output_filename, calc, plot)
            .depends_on(self.load())
            .use_cores(-1)
        )

    def register_qc_splicing(self):
        """How many reads were spliced? How many of those splices were known splice sites,
        how many were novel"""
        output_filename = self.result_dir / f"{self.name}_splice_sites.png"

        def calc():
            from mbf_bam import count_introns

            bam_filename, bam_index_name = self.get_bam_names()
            counts_per_chromosome = count_introns(bam_filename, bam_index_name)
            known_splice_sites_by_chr = {
                chr: set() for chr in self.genome.get_chromosome_lengths()
            }
            for gene in self.genome.genes.values():
                for start, stop in zip(*gene.introns_all):
                    known_splice_sites_by_chr[gene.chr].add((start, stop))
            total_counts = collections.Counter()
            known_count = 0
            unknown_count = 0
            for chr, counts in counts_per_chromosome.items():
                for k, v in counts.items():
                    if k[0] == 0xFFFFFFFF:
                        intron_counts = 0xFFFFFFFF - k[1]
                        total_counts[intron_counts] += v
                    else:
                        if k in known_splice_sites_by_chr[chr]:
                            known_count += v
                        else:
                            unknown_count += v
            result = {"side": [], "x": [], "count": []}
            result["side"].append("splice sites")
            result["x"].append("unknown")
            result["count"].append(unknown_count)
            result["side"].append("splice sites")
            result["x"].append("known")
            result["count"].append(known_count)

            for x, count in total_counts.items():
                result["side"].append("reads with x splices")
                result["x"].append(x)
                result["count"].append(count)

            return pd.DataFrame(result)

        def plot(df):
            return (
                dp(df)
                .p9()
                .theme_bw()
                .add_bar("x", "count", stat="identity")
                .facet_wrap("side", scales="free", ncol=1)
                .scale_y_continuous(labels=lambda xs: ["%.2g" % x for x in xs])
                .title(self.name)
                .theme(panel_spacing_y=0.2)
                .render(output_filename)
            )

        return register_qc(
            ppg.PlotJob(output_filename, calc, plot)
            .depends_on(self.load())
            .use_cores(-1)
        )


__all__ = [Sample, AlignedSample]
