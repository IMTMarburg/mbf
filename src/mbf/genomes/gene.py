import attr
import numpy as np
from mbf_nested_intervals import IntervalSet
from .common import reverse_complement


@attr.s(slots=True)
class Gene:
    gene_stable_id = attr.ib()
    name = attr.ib()
    chr = attr.ib()
    start = attr.ib()
    stop = attr.ib()
    strand = attr.ib()
    biotype = attr.ib()
    transcripts = attr.ib()
    genome = attr.ib()

    @property
    def tss(self):
        return self.start if self.strand == 1 else self.stop

    @property
    def tes(self):
        return self.start if self.strand != 1 else self.stop

    @property
    def introns_strict(self):
        """Get truly intronic regions - ie. not covered by any exon for this gene
        result is a a tuple of np arrays, (starts, stops)
        By it's definition, the introns are disjunct
        """
        gene_start = self.start
        gene_stop = self.stop
        exons = []
        for tr in self.transcripts:
            try:
                exons.extend(tr.exons)
            except TypeError:  # pragma: no cover
                raise ValueError(f"No exons defined for {tr.transcript_stable_id}")
        return IntervalSet.from_tuples(exons).invert(gene_start, gene_stop).to_numpy()

    @property
    def introns_all(self):
        """Get intronic regions - ie. an intron in any of the transcripts.
        May contain repetitions and overlaps and is not sorted!
        """
        gene_start = self.start
        gene_stop = self.stop
        introns = [], []
        for tr in self.transcripts:
            try:
                starts, stops = (
                    IntervalSet.from_tuples(tr.exons)
                    .invert(gene_start, gene_stop)
                    .to_numpy()
                )
            except TypeError:  # pragma: no cover
                raise ValueError(f"No exons defined for {tr.transcript_stable_id}")
            introns[0].extend(starts)
            introns[1].extend(stops)
        return introns

    @property
    def _exons(self):
        """Common code to exons_merged and exons_overlapping"""
        exons = []
        for tr in self.transcripts:
            try:
                exons.extend(tr.exons)
            except TypeError:  # pragma: no cover
                raise ValueError(f"No exons defined for {tr.transcript_stable_id}")
        return exons

    @property
    def exons_merged(self):
        """Get the merged exon regions for a gene given by gene_stable_id
        result is a a tuple of np arrays, (starts, stops)
        """
        return IntervalSet.from_tuples(self._exons).merge_hull().to_numpy()

    @property
    def exons_overlapping(self):
        """Get the overlapping exon regions for a gene given by gene_stable_id
        result is a a tuple of np arrays, (starts, stops)
        not sorted
        """
        return self._reformat_exons(self._exons)

    def _reformat_exons(self, exons):
        """Turn exons [(start, stop), ...] into [[start, ...], [stop, ...]
        """
        exons.sort()
        return np.array([x[0] for x in exons]), np.array([x[1] for x in exons])

    @property
    def _exons_protein_coding(self):
        """common code for the exons_protein_coding_* propertys"""
        exons = []
        for tr in self.transcripts:
            if tr.biotype == "protein_coding":
                exons.extend(tr.exons)
        return exons

    @property
    def exons_protein_coding_merged(self):
        """Get the merged exon regions for a gene , only for protein coding exons.
        Empty result on non protein coding genes
        result is a a tuple of np arrays, (starts, stops)
        """
        return (
            IntervalSet.from_tuples(self._exons_protein_coding).merge_hull().to_numpy()
        )

    @property
    def exons_protein_coding_overlapping(self):
        """Get the overlapping exon regions for a gene, only for protein coding transcripts.
        Empty result on non protein coding genes

        Result is a DataFrame{chr, strand, start, stop}

        We test biotype on transcripts, not on genes,
        because for example polymorphismic_pseudogenes can have protein coding variants.
        """
        return self._reformat_exons(self._exons_protein_coding)

    def get_reads_in_exons(self, bam, dedup="gene"):
        """Return a list of all the reads in the bam that
        have an overlap with the exon regions

        dedup may be:
            'gene' - a given read may only count once in this gene (default)
            False - no dedup (but each alignment only once even if it spans multiple exons)
            'primary_only': - only return reads with tag HI==1
        """
        result = {}
        if dedup == "gene":
            for start, stop in zip(*self.exons_merged):
                for r in bam.fetch(self.chr, start, stop):
                    if r.get_overlap(start, stop) > 0:
                        result[r.query_name] = r
        elif dedup is False:
            for start, stop in zip(*self.exons_merged):
                for r in bam.fetch(self.chr, start, stop):
                    if r.get_overlap(start, stop) > 0:
                        # which amounts to no-dedud except if the aligner is
                        # *really* buggy
                        result[len(result)] = r
        elif dedup == "primary_only":
            for start, stop in zip(*self.exons_merged):
                for r in bam.fetch(self.chr, start, stop):
                    if r.get_overlap(start, stop) > 0 and r.get_tag("HI") == 1:
                        # which amounts to no-dedud except if the aligner is
                        # *really* buggy
                        result[r.query_name, r.pos, r.cigarstring] = r
        else:
            raise ValueError("Invalid dedup value")

        return list(result.values())


@attr.s(slots=True)
class Transcript:
    transcript_stable_id = attr.ib()
    gene_stable_id = attr.ib()
    name = attr.ib()
    chr = attr.ib()
    start = attr.ib()
    stop = attr.ib()
    strand = attr.ib()
    biotype = attr.ib()
    exons = attr.ib()
    exon_stable_ids = attr.ib()
    gene = attr.ib()
    genome = attr.ib()

    @property
    def exons_tuples(self):
        return [(start, stop) for (start, stop) in self.exons]

    @property
    def introns(self):
        """Return [(start, stop),...] for all introns in the transcript
        Order is in genomic order.
        Intron is defined as everything inside tss..tes that is not an exon,
        so if a gene, by any reason would extend beyond it's exons,
        that region would also be covered.
        """
        gene_start = self.gene.start
        gene_stop = self.gene.stop
        exons = sorted(self.exons_tuples)
        return IntervalSet.from_tuples(exons).invert(gene_start, gene_stop).to_tuples()

    @property
    def cdna(self):
        """Get the cdna sequence as defined by cdna.fasta"""
        return self.genome.get_cdna_sequence(self.transcript_stable_id)

    @property
    def mrna(self):
        """The mRNA sequence after splicing.
        (forward strand - ie. ATG is ATG)

        unlike cdna, this is build dynamically from the genome_sequence
        and exon definition and is available for non-protein coding transcripts
        """
        seq = "".join(
            [
                self.genome.get_genome_sequence(self.chr, start, stop)
                for (start, stop) in self.exons
            ]
        )
        if self.strand == -1:
            seq = reverse_complement(seq)
        return seq

    @property
    def coordinate_translations(self):
        """Return a [genomic_coordinates]
        the transcript coordinate relative to the  5' start
        is the index
        """
        result = []
        tr_pos = 0
        if self.strand == 1:
            for exon_start, exon_stop in self.exons:
                for ii in range(exon_start, exon_stop):
                    # result.append((tr_pos, ii))
                    result.append(ii)
                    tr_pos += 1
        else:
            for exon_start, exon_stop in reversed(self.exons):
                for ii in range(exon_stop, exon_start, -1):
                    # result.append((tr_pos, ii - 1))
                    result.append(ii - 1)
                    tr_pos += 1
        return result
