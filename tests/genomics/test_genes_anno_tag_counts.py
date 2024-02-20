import pytest
from pathlib import Path
import numpy as np
import pandas as pd
from mbf.genomics.genes.anno_tag_counts import IntervalStrategyIntron

from mbf.genomics.genes import Genes, anno_tag_counts
from mbf.qualitycontrol import prune_qc, get_qc_jobs, qc_disabled
from mbf.qualitycontrol.testing import assert_image_equal
from .shared import MockGenome, force_load, run_pipegraph, RaisesDirectOrInsidePipegraph
from .old_reference_code import NormalizationCPMBiotypes, NormalizationTPMBiotypes
import pysam


def MockBam(chr, length, mode=1):
    """a number of reads of length 1 at any position!"""
    chrs = [chr]
    fn = "mock.bam"
    bam = pysam.Samfile(
        fn, "wb", reference_names=chrs, reference_lengths=[100_000] * len(chrs)
    )

    def write_read(pos, strand, name, nh=1):
        al = pysam.AlignedSegment()
        al.pos = pos
        al.reference_id = 0
        al.is_reverse = strand == -1
        al.query_name = name
        al.seq = "A"
        al.cigar = ((0, 1),)
        al.tags = (("NH", nh),)
        bam.write(al)

    for ii in range(0, length):
        if mode == 0:
            write_read(ii, 1, "fw_%i" % ii)
        elif mode == 1:
            write_read(ii, 1, "fw_%i" % ii)
            write_read(ii, -1, "rv_%i" % ii)
        elif mode == 2:
            write_read(ii, 1, "fw_%i" % ii)
            write_read(ii, -1, "rv_%i" % ii)
            write_read(ii, -1, "rv2_%i" % ii)
        elif mode == 3:  # good for testing the 'count one read in one gene once'
            write_read(ii, 1, "fw_%i" % ii, 4)
            write_read(ii, 1, "fw_%i" % ii, 4)
            write_read(ii, 1, "fw_%i" % ii, 4)
            write_read(ii, 1, "fw_%i" % ii, 4)
            write_read(ii, -1, "rv2_%i" % ii)
        elif mode == 4:
            write_read(ii, 1, "fw_%i_a" % ii)
            write_read(ii, 1, "fw_%i_b" % ii)
            write_read(ii, 1, "fw_%i_c" % ii)
            write_read(ii, 1, "fw_%i_d" % ii)
            write_read(ii, -1, "rv2_%i_a" % ii)
        else:
            raise ValueError("Invalid mode")
    bam.close()
    pysam.sort("-o", "mock_sorted.bam", fn)
    pysam.index("mock_sorted.bam")
    res = pysam.Samfile("mock_sorted.bam", "rb")
    return res


def MockBamFixed(reads):
    """Each read is a
    dict: strand, regions, tags, qname, is_read1
    """
    ii = 0
    for read in reads:
        if "qname" not in read:
            name = "read_%i" % ii
            ii = +1
            read["qname"] = name
        if "chr" not in read:
            read["chr"] = "1"
        if "tags" not in read:
            read["tags"] = {}
        if "is_read1" not in read:
            read["is_read1"] = True
    chrs = sorted(set([r["chr"] for r in reads]))
    fn = "mock.bam"
    if Path(fn).exists():
        Path(fn).unlink()
    bam = pysam.Samfile(
        fn, "wb", reference_names=chrs, reference_lengths=[100_000] * len(chrs)
    )
    for read in reads:
        al = pysam.AlignedSegment()
        al.query_name = read["qname"]
        read_len = sum([x[1] - x[0] for x in read["regions"]])
        al.query_sequence = "A" * read_len
        al.reference_id = bam.references.index(read["chr"])
        al.is_read1 = read["is_read1"]
        al.is_reverse = read["strand"] == -1
        t = [tuple(x) for x in read["tags"].items()]
        al.tags = t
        al.pos = read["regions"][0][0]
        last_stop = None
        cigar = []
        for start, stop in read["regions"]:
            if last_stop:
                cigar.append((2, start - last_stop))
            last_stop = stop
            cigar.append((0, stop - start))
        al.cigar = tuple(cigar)
        bam.write(al)
    bam.close()
    pysam.sort("-o", "mock_sorted.bam", fn)
    pysam.index("mock_sorted.bam")
    res = pysam.Samfile("mock_sorted.bam", "rb")
    return res


def FakePaireEndBam(chr, length, reverse=False):
    chrs = [chr]
    fn = "mock.bam"
    bam = pysam.Samfile(
        fn, "wb", reference_names=chrs, reference_lengths=[100_000] * len(chrs)
    )

    def write_read(pos, strand, name):
        al = pysam.AlignedSegment()
        al.pos = pos
        al.reference_id = 0
        al.is_reverse = strand == -1
        al.query_name = name
        al.seq = "A"
        al.cigar = ((0, 1),)
        al.tags = (("NH", 1),)
        bam.write(al)

    for ii in range(0, length):
        if not reverse:
            write_read(ii, 1, "read_%i" % ii, is_read1=True)
            write_read(ii, -1, "read_%i" % ii, is_read1=False)
        else:
            write_read(ii, 1, "read_%i" % ii, is_read1=False)
            write_read(ii, -1, "read_%i" % ii, is_read1=True)

    bam.close()
    pysam.sort("-o", "mock_sorted.bam", fn)
    pysam.index("mock_sorted.bam")
    res = pysam.Samfile("mock_sorted.bam", "rb")
    return res


class MockLane(object):
    def __init__(self, bam, genome, name="fakelane1"):
        self.name = name
        self.genome = genome
        self.bam_filename = bam.filename
        self.bam_index_filename = bam.index_filename
        self.vid = "FAkeVID"

    def get_bam(self):
        return pysam.Samfile(self.bam_filename, self.bam_filename_index_filename)

    def get_bam_names(self):
        bam_filename = self.bam_filename.decode("utf-8")
        bam_index_name = self.bam_index_filename
        if bam_index_name is None:
            bam_index_name = bam_filename + ".bai"
        else:
            bam_index_name = str(bam_index_name)

        return bam_filename, bam_index_name

    def load(self):
        pass


@pytest.mark.usefixtures("both_ppg_and_no_ppg_no_qc")
class TestGeneCount:
    def test_unstranded(self):

        bam = MockBam("1", 10000)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 1000,
                    "tes": 2000,
                    "chr": "1",
                }
            ]
        )
        transcripts = pd.DataFrame({"name": [], "chr": [], "exons": [], "strand": []})
        genome = MockGenome(genes, transcripts, {"1": 10000})
        genes = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.GeneUnstranded(lane)
        genes += anno
        force_load(genes.annotate)
        run_pipegraph()
        df = genes.df
        assert (df[anno.columns[0]] == np.array([2000])).all()

    def test_stranded_raises_on_too_many_reverse_reads(self):

        bam = MockBam("1", 10000, mode=2)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 1000,
                    "tes": 2000,
                    "chr": "1",
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                },
            ]
        )
        transcripts = pd.DataFrame({"name": [], "chr": [], "exons": [], "strand": []})
        genome = MockGenome(genes, transcripts, {"1": 10000})
        genes = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.GeneStranded(lane)
        with RaisesDirectOrInsidePipegraph(ValueError):
            force_load(genes.add_annotator(anno))

    def test_stranded(self):

        bam = MockBam("1", 10000, mode=2)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 2000,
                    "tes": 1000,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeC",
                    "name": "C",
                    "tss": 400,
                    "tes": 410,
                    "chr": "1",
                    "strand": 1,
                },
                # to get around the too many reverse reads threshold
                {
                    "gene_stable_id": "FakeD",
                    "name": "D",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeD2",
                    "name": "D2",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeD3",
                    "name": "D3",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeD4",
                    "name": "D4",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeD5",
                    "name": "D5",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeD6",
                    "name": "D6",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
                {
                    "gene_stable_id": "FakeD7",
                    "name": "D7",
                    "tss": 2010,
                    "tes": 1910,
                    "chr": "1",
                    "strand": -1,
                },
            ]
        ).sort_values("gene_stable_id")

        transcripts = pd.DataFrame({"name": [], "chr": [], "exons": [], "strand": []})
        genome = MockGenome(genes, transcripts, {"1": 10000})
        genes = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.GeneStranded(lane)
        force_load(genes.add_annotator(anno))
        run_pipegraph()
        assert list(genes.df.index) == list(
            range(len(genes.df))
        ), "make sure the Genes.df has an index 0..n"  # make sure that the Genes have an index that goes 0..n

        assert (
            genes.df.gene_stable_id
            == [
                "FakeC",
                "FakeA",
                "FakeB",
                "FakeD",
                "FakeD2",
                "FakeD3",
                "FakeD4",
                "FakeD5",
                "FakeD6",
                "FakeD7",
            ]
        ).all()

        assert (
            genes.df[anno.columns[0]]
            == np.array([10, 2000, 200, 200, 200, 200, 200, 200, 200, 200])
        ).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg_no_qc")
class TestExonSmartCount:
    def test_validation(self):
        lane = object()
        with pytest.raises(ValueError):
            anno_tag_counts.ExonSmartUnstranded(lane)

    def test_unstranded(self):

        bam = MockBam("1", 10000, mode=2)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 500,
                    "tes": 1000,
                    "chr": "1",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "transcript_stable_id": "fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "transcript_stable_id": "fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "transcript_stable_id": "fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "transcript_stable_id": "fakeab3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartUnstranded(lane)
        genes = Genes(genome)
        force_load(genes.add_annotator(anno))
        run_pipegraph()
        df = genes.df.sort_values("name")
        print(df)
        print(genome.genes["FakeA"].exons_protein_coding_overlapping)
        assert (df[anno.columns[0]] == np.array([300, (100 + 150) * 3])).all()

    def test_stranded(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 500,
                    "tes": 600,
                    "chr": "1",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane, "output.bam")
        anno.count_strategy.disable_sanity_check = True
        force_load(gr.add_annotator(anno))
        run_pipegraph()
        df = gr.df.sort_values("name")
        print(df)
        assert (df[anno.columns[0]] == np.array([100 * 4, (100 + 150) * 1])).all()
        import pysam

        f = pysam.Samfile("output.bam")
        reads = list(f.fetch())
        assert (
            len(reads) == gr.df[anno.columns[0]].sum()
        )  # since we have each fw read 4 times in the bam, and we keep all of them

    def test_stranded_dedup(self):

        bam = MockBam("1", 10000, mode=3)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 500,
                    "tes": 700,
                    "chr": "1",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane)
        force_load(gr.add_annotator(anno))
        run_pipegraph()
        assert (gr.df[anno.columns[0]] == np.array([100 * 1, (100 + 150) * 1])).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg_no_qc")
class TestExonCount:
    def test_unstranded(self):

        bam = MockBam("1", 10000, mode=2)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 500,
                    "tes": 2000,
                    "chr": "1",
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonUnstranded(lane)
        force_load(gr.add_annotator(anno))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")
        assert (df[anno.columns[0]] == np.array([300, (1100) * 3])).all()

    def test_stranded(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonStranded(lane)
        anno.count_strategy.disable_sanity_check = True
        force_load(gr.add_annotator(anno))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")
        assert (df[anno.columns[0]] == np.array([100 * 4, (1100) * 1])).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg_no_qc")
class TestCPM:
    def test_stranded(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane)
        anno.count_strategy.disable_sanity_check = True
        anno2 = anno_tag_counts.NormalizationCPM(anno)
        force_load(gr.add_annotator(anno2))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")
        assert (df[anno.columns[0]] == np.array([100 * 4, (250) * 1])).all()
        assert (
            df[anno2.columns[0]] == np.array([400 * 1e6 / 650.0, 250 * 1e6 / 650.0])
        ).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg_no_qc")
class TestTPM:
    def test_stranded(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                    "strand": 1,
                },  # length: 1000
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "strand": -1,
                },  # length 2000
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                # length 100
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },  # length = 200
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },  # 1100...
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },  # length: 250
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane)
        anno.count_strategy.disable_sanity_check = True
        anno2 = anno_tag_counts.NormalizationTPM(anno)

        force_load(gr.add_annotator(anno2))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")
        tpm_factor = 1e6 / (400.0 / 100 + 250.0 / 250)

        assert (
            df[anno2.columns[0]]
            == np.array([400 / 100.0 * tpm_factor, 250 / 250.0 * tpm_factor])
        ).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg_no_qc")
class TestCPMBiotypes:
    def test_stranded(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "protein_coding",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "misc_RNA",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane)
        anno.count_strategy.disable_sanity_check = True
        anno2 = NormalizationCPMBiotypes(anno, ("protein_coding",))
        force_load(gr.add_annotator(anno2))
        gr_biotype = gr.filter("pc", lambda df: df["biotype"] == "protein_coding")
        anno3 = anno_tag_counts.NormalizationCPM(anno)
        force_load(gr_biotype.add_annotator(anno3))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")

        assert df[anno2.columns[0]][0] == 1e6
        assert np.isnan(df[anno2.columns[0]][1])

        assert gr_biotype.df[anno3.columns[0]][0] == 1e6
        assert len(gr_biotype.df) == 1

    def test_stranded2(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "protein_coding",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "misc_RNA",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane)
        anno.count_strategy.disable_sanity_check = True
        anno2 = NormalizationCPMBiotypes(anno, ("protein_coding", "misc_RNA"))
        force_load(gr.add_annotator(anno2))
        gr_biotype = gr.filter(
            "pc", lambda df: df["biotype"].isin(["protein_coding", "misc_RNA"])
        )
        anno3 = anno_tag_counts.NormalizationCPM(anno)
        force_load(gr_biotype.add_annotator(anno3))
        force_load(gr_biotype.add_annotator(anno2))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")

        assert (
            df[anno2.columns[0]] == np.array([400 * 1e6 / 650, 250 * 1e6 / 650])
        ).all()
        assert (
            gr_biotype.df[anno3.columns[0]] == gr_biotype.df[anno2.columns[0]]
        ).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg_no_qc")
class TestTPMBiotypes:
    def test_stranded(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "protein_coding",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "misc_RNA",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane)
        anno.count_strategy.disable_sanity_check = True
        anno2 = NormalizationTPMBiotypes(anno, ("protein_coding",))
        force_load(gr.add_annotator(anno2))
        gr_biotype = gr.filter("pc", lambda df: df["biotype"].isin(["protein_coding"]))
        anno3 = anno_tag_counts.NormalizationTPM(anno)
        force_load(gr_biotype.add_annotator(anno3))
        force_load(gr_biotype.add_annotator(anno2))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")
        # length_by_gene = anno.interval_strategy.get_interval_lengths_by_gene(genome)
        # counts = np.array([400, 250.0])
        # lengths = [100.0, 250]
        # xi = counts / lengths
        # tpms = xi * (1e6 / xi.sum())
        assert df[anno2.columns[0]][0] == 1e6
        assert np.isnan(df[anno2.columns[0]][1])
        assert (
            gr_biotype.df[anno3.columns[0]][0] == gr_biotype.df[anno2.columns[0]][0]
        ).all()

    def test_stranded2(self):

        bam = MockBam("1", 10000, mode=4)
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "protein_coding",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "misc_RNA",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        gr = Genes(genome)
        lane = MockLane(bam, genome)
        anno = anno_tag_counts.ExonSmartStranded(lane)
        anno.count_strategy.disable_sanity_check = True
        anno2 = NormalizationTPMBiotypes(anno, ("protein_coding", "misc_RNA"))
        force_load(gr.add_annotator(anno2))
        gr_biotype = gr.filter(
            "pc", lambda df: df["biotype"].isin(["protein_coding", "misc_RNA"])
        )
        anno3 = anno_tag_counts.NormalizationTPM(anno)
        force_load(gr_biotype.add_annotator(anno3))
        force_load(gr_biotype.add_annotator(anno2))
        run_pipegraph()
        df = gr.df.sort_values("gene_stable_id")

        length_by_gene = anno.interval_strategy.get_interval_lengths_by_gene(genome)
        tpm_factor = 1e6 / (400.0 / 100 + 250.0 / 250)
        print(length_by_gene)
        assert (
            df[anno2.columns[0]]
            == np.array([400 / 100.0 * tpm_factor, 250 / 250.0 * tpm_factor])
        ).all()
        assert (
            gr_biotype.df[anno3.columns[0]][0] == gr_biotype.df[anno2.columns[0]][0]
        ).all()

    def test_interval_strategy_intron(self):
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 0,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "protein_coding",
                    "strand": 1,
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 4000,
                    "tes": 2000,
                    "chr": "1",
                    "biotype": "misc_RNA",
                    "strand": -1,
                },
            ]
        )
        transcripts = pd.DataFrame(
            [
                {
                    "name": "Fakea1",
                    "chr": "1",
                    "exons": [(500, 600)],
                    "strand": 1,
                    "gene_stable_id": "FakeA",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb1",
                    "chr": "1",
                    "exons": [(3900, 4000), (2900, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
                {
                    "name": "Fakeb2",
                    "chr": "1",
                    "exons": [(3900, 4000), (2000, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "misc_RNA",
                },
                {
                    "name": "Fakeb3",
                    "chr": "1",
                    "exons": [(3900, 4000), (2850, 3000)],
                    "strand": -1,
                    "gene_stable_id": "FakeB",
                    "biotype": "protein_coding",
                },
            ]
        )
        genome = MockGenome(genes, transcripts, {"1": 10000})
        i = IntervalStrategyIntron()
        actual = i._get_interval_tuples_by_chr(genome)
        assert "1" in actual
        assert len(actual["1"]) == 2
        assert len(actual["1"][0][2]) == 2


class TestFPKM:
    def test_fpkm_raises(self):
        with pytest.raises(NotImplementedError):
            anno_tag_counts.NormalizationFPKM(55)


@pytest.mark.usefixtures("new_pipegraph_no_qc")
class TestPPG:
    def test_cores_needed(self):
        genes = pd.DataFrame(
            [
                {
                    "gene_stable_id": "FakeA",
                    "name": "A",
                    "tss": 15,
                    "tes": 300,
                    "chr": "1",
                },
                {
                    "gene_stable_id": "FakeB",
                    "name": "B",
                    "tss": 1200,
                    "tes": 1000,
                    "chr": "1",
                },
            ]
        )
        import pypipegraph as ppg

        if not ppg.inside_ppg():
            raise ValueError()
        transcripts = pd.DataFrame({"name": [], "chr": [], "exons": [], "strand": []})
        genome = MockGenome(genes, transcripts, {"1": 10000})
        bam = MockBam("1", 100)
        lane = MockLane(bam, genome, name="shu")
        g = Genes(genome)
        annos = [
            anno_tag_counts.ExonSmartStranded(lane),
            anno_tag_counts.ExonSmartUnstranded(lane),
        ]
        for anno in annos:
            g.add_annotator(anno)
        # ppg.run_pipegraph()
        g.load_strategy.fix_anno_tree()

        for anno in annos:
            j = g.anno_jobs[anno.get_cache_name()]
            assert j.lfg.cores_needed == -1

    def test_tpm_subset_of_genes(self, new_pipegraph_no_qc):
        from mbf.sampledata import get_human_22_fake_genome, get_sample_data
        from mbf.qualitycontrol import disable_qc
        import pypipegraph as ppg
        import mbf.align
        from mbf.genomics.annotator import Annotator

        ppg.util.global_pipegraph.quiet = False

        genome = get_human_22_fake_genome()
        lane = mbf.align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        genes = Genes(genome)
        pc = genes.filter("protein_coding", lambda df: df.biotype == "protein_coding")
        raw = anno_tag_counts.ExonSmartStranded(lane)
        tpm = anno_tag_counts.NormalizationTPM(raw)
        pc.add_annotator(tpm)
        pc.write()
        ppg.run_pipegraph()
        assert pc.df[tpm.columns[0]].sum() == pytest.approx(1e6)
        assert not raw.columns[0] in genes.df.columns
        with_filtered = pc.df[tpm.columns[0]]

        new_pipegraph_no_qc.new_pipegraph()
        disable_qc()
        genome = get_human_22_fake_genome()
        lane = mbf.align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic

        genes = Genes(genome)
        pc = genes.filter("protein_coding2", lambda df: df.biotype == "protein_coding")
        raw = anno_tag_counts.ExonSmartStranded(lane)

        class RawToNanOnBiotype(Annotator):
            def __init__(self):
                self.columns = ["Rawww"]
                self.genome = genome
                self.interval_strategy = raw.interval_strategy

            def calc(self, df):
                column = df[raw.columns[0]].copy()
                column[df["biotype"] != "protein_coding"] = np.nan
                return column

            def dep_annos(self):
                return [raw]

        tpm = anno_tag_counts.NormalizationTPM(RawToNanOnBiotype())
        genes.add_annotator(tpm)
        pc.write()
        ppg.run_pipegraph()
        assert (pc.df[tpm.columns[0]] == with_filtered).all()


@pytest.mark.usefixtures("new_pipegraph")
class TestQC:
    def test_qc_distribution(self):
        from mbf.sampledata import get_human_22_fake_genome, get_sample_data
        import mbf.align

        genome = get_human_22_fake_genome()
        lane = mbf.align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        lane2 = mbf.align.AlignedSample(
            "test_lane2",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        genes = Genes(genome)
        for lane in [lane, lane2]:
            anno = anno_tag_counts.GeneStranded(lane)
            genes += anno
            genes += anno_tag_counts.NormalizationTPM(anno)
            genes += anno_tag_counts.NormalizationCPM(anno)
        assert not qc_disabled()
        prune_qc(lambda job: "read_distribution" in job.job_id)
        run_pipegraph()
        qc_jobs = list(get_qc_jobs())
        qc_jobs = [x for x in qc_jobs if not x._pruned]
        assert len(qc_jobs) == 4  # three from our annos, one gene_unstranded for
        cpm_job = [x for x in qc_jobs if "CPM" in str(x.filenames[0])][0]
        tpm_job = [x for x in qc_jobs if "TPM" in str(x.filenames[0])][0]
        raw_job = [
            x
            for x in qc_jobs
            if (
                not "TPM" in str(x.filenames[0])
                and not "CPM" in str(x.filenames[0])
                and not "unstranded" in str(x.filenames[0])
            )
        ][0]
        assert_image_equal(raw_job.filenames[0])
        assert_image_equal(tpm_job.filenames[0], "_tpm")
        assert_image_equal(cpm_job.filenames[0], "_cpm")

    def test_qc_distribution_single_gene_sequenced(self):
        from mbf.sampledata import get_human_22_fake_genome, get_sample_data
        import mbf.align
        import pypipegraph as ppg

        genome = get_human_22_fake_genome()
        lane = mbf.align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        print("after creation freeze", ppg.util.freeze(lane))
        lane2 = mbf.align.AlignedSample(
            "test_lane2",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        genes = Genes(genome)

        def fake_calc(df):
            result = np.zeros((len(df),), float)
            result[2] = 1
            return pd.Series(result)

        for lane in [lane, lane2]:
            print("adding", lane.name)
            print(genes.annotators)
            print("just before freeze", ppg.util.freeze(lane))
            anno = anno_tag_counts.GeneUnstranded(lane)
            anno.calc = fake_calc
            genes += anno
            print("done adde", lane.name)

        assert not qc_disabled()
        prune_qc(lambda job: "read_distribution" in job.job_id)
        run_pipegraph()
        qc_jobs = list(get_qc_jobs())
        qc_jobs = [x for x in qc_jobs if not x._pruned]
        assert len(qc_jobs) == 1  # three from our annos, one gene_unsrtanded for

        raw_job = qc_jobs[0]
        assert_image_equal(raw_job.filenames[0])

    def test_qc_distribution_no_data(self):
        from mbf.sampledata import get_human_22_fake_genome, get_sample_data
        import mbf.align

        genome = get_human_22_fake_genome()
        lane = mbf.align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        lane2 = mbf.align.AlignedSample(
            "test_lane2",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        genes = Genes(genome)
        genes = genes.filter("none", lambda df: [False] * len(df))

        def fake_calc(df):
            result = np.zeros((len(df),), float)
            return pd.Series(result)

        for lane in [lane, lane2]:
            anno = anno_tag_counts.GeneUnstranded(lane)
            anno.calc = fake_calc
            genes += anno

        assert not qc_disabled()
        prune_qc(lambda job: "read_distribution" in job.job_id and "none" in job.job_id)
        run_pipegraph()
        qc_jobs = list(get_qc_jobs())
        qc_jobs = [x for x in qc_jobs if not x._pruned]
        assert len(qc_jobs) == 1  # three from our annos, one gene_unsrtanded for

        raw_job = qc_jobs[0]
        assert_image_equal(raw_job.filenames[0])

    def test_qc_pca(self):
        import mbf.sampledata

        ddf, a, b = mbf.sampledata.get_pasilla_data_subset()
        annos = []
        for x in a + b:
            anno = anno_tag_counts.NormalizationCPM(x)
            ddf += anno
            annos.append(anno)
        ddf.write()
        prune_qc(lambda job: "pca" in job.job_id)
        run_pipegraph()
        qc_jobs = list(get_qc_jobs())
        qc_jobs = [x for x in qc_jobs if not x._pruned]
        assert len(qc_jobs) == 1
        assert_image_equal(qc_jobs[0].filenames[0])

    def test_qc_pca_filtered(self):
        import mbf.sampledata

        ddf, a, b = mbf.sampledata.get_pasilla_data_subset()
        annos = []
        for x in a + b:
            anno = anno_tag_counts.NormalizationCPM(x)
            ddf += anno
            annos.append(anno)
        ddf2 = ddf.filter("filtered", lambda df: df[anno.columns[0]] >= 100, [anno])
        ddf2.write()
        prune_qc(lambda job: "pca" in job.job_id)
        run_pipegraph()
        qc_jobs = list(get_qc_jobs())
        qc_jobs = [x for x in qc_jobs if not x._pruned]
        assert len(qc_jobs) == 2
        assert_image_equal(qc_jobs[0].filenames[0])
        assert_image_equal(qc_jobs[1].filenames[0], "_filtered")

    def test_qc_pca_single_sample(self):
        import mbf.sampledata

        ddf, a, b = mbf.sampledata.get_pasilla_data_subset()
        annos = []
        for x in a + b:
            anno = anno_tag_counts.NormalizationCPM(x)
            ddf += anno
            annos.append(anno)
            break
        prune_qc(lambda job: "pca" in job.job_id)
        ddf.write()
        run_pipegraph()
        qc_jobs = list(get_qc_jobs())
        qc_jobs = [x for x in qc_jobs if not x._pruned]
        assert len(qc_jobs) == 1
        assert_image_equal(qc_jobs[0].filenames[0])

    def test_qc_pca_no_data(self):
        import mbf.sampledata

        ddf, a, b = mbf.sampledata.get_pasilla_data_subset()
        ddf2 = ddf.filter("no_data", lambda df: [False] * len(df))
        annos = []
        for x in a + b:
            anno = anno_tag_counts.NormalizationCPM(x)
            ddf2 += anno
            annos.append(anno)
        ddf2.write()
        prune_qc(lambda job: "pca" in job.job_id)
        run_pipegraph()
        qc_jobs = list(get_qc_jobs())
        qc_jobs = [x for x in qc_jobs if not x._pruned]
        assert len(qc_jobs) == 1
        assert_image_equal(qc_jobs[0].filenames[0])
