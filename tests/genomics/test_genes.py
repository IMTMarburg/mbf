import pypipegraph as ppg
import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from pandas.testing import assert_frame_equal
from mbf.fileformats.bed import read_bed
import mbf.genomics.regions as regions
import mbf.genomics.genes as genes
from mbf.genomics.annotator import Constant

from .shared import (
    get_genome,
    get_genome_chr_length,
    force_load,
    run_pipegraph,
    RaisesDirectOrInsidePipegraph,
    MockGenome,
)


@pytest.mark.usefixtures("new_pipegraph")
class TestGenesLoadingPPGOnly:
    def test_loading_from_genome_is_singletonic(self):
        genome = get_genome()
        print(genome)
        genesA = genes.Genes(genome)
        genesB = genes.Genes(genome)
        assert genesA is genesB
        filterA = genesA.filter("fa", lambda df: df.index[:10])
        filterAa = genesA.filter("faa", lambda df: df.index[:10])
        filterB = genesB.filter("fab", lambda df: df.index[:10])
        assert not (filterA is genesA)
        assert not (filterAa is filterA)
        assert not (filterAa is filterB)
        with pytest.raises(ValueError):  # can't have a different loading func
            filterB = genesB.filter("fab", lambda df: df.index[:15])
        force_load(filterA.load)
        ppg.run_pipegraph()
        assert len(filterA.df) == 10


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenesLoading:
    def test_basic_loading_from_genome(self):
        g = genes.Genes(get_genome())
        force_load(g.load())
        run_pipegraph()
        assert len(g.df) == 246
        assert (g.df["gene_stable_id"][:3] == ["CRP_001", "CRP_002", "CRP_003"]).all()
        assert g.df["gene_stable_id"].iloc[-1] == "CRP_182"
        assert g.df["start"].iloc[-1] == 158_649 - 1
        assert g.df["stop"].iloc[-1] == 159_662
        assert g.df["strand"].iloc[-1] == -1

    def test_filtering_with_annotator(self):
        import mbf.genomics

        g = genes.Genes(get_genome())

        class CopyAnno(mbf.genomics.annotator.Annotator):
            def __init__(self):
                self.columns = ["copy"]

            def calc(self, df):
                return pd.DataFrame({"copy": df["gene_stable_id"]})

        g += CopyAnno()
        filtered = g.filter("a", ("gene_stable_id", "==", "CRP_003"))
        force_load(filtered.annotate())
        run_pipegraph()
        assert (filtered.df["gene_stable_id"] == ["CRP_003"]).all()
        assert (filtered.df["copy"] == ["CRP_003"]).all()

    def test_alternative_loading_raises_on_non_df(self):
        with RaisesDirectOrInsidePipegraph(ValueError):
            g = genes.Genes(get_genome_chr_length(), lambda: None, "myname")
            force_load(g.load())

    def test_alternative_loading_raises_on_missing_column(self, both_ppg_and_no_ppg):
        df = pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 5000,
                    "tes": 5500,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        )

        def inner_tss():
            df2 = df.copy()
            df2 = df2.drop("tss", axis=1)
            g = genes.Genes(get_genome(), lambda: df2, name="sha")
            g.load()
            # run_pipegraph()

        def inner_chr():
            df2 = df.copy()
            df2 = df2.drop("chr", axis=1)
            g = genes.Genes(get_genome(), lambda: df2, name="shu")
            g.load()
            # run_pipegraph()

        def inner_tes():
            df2 = df.copy()
            df2 = df2.drop("tes", axis=1)
            g = genes.Genes(get_genome(), lambda: df2, name="shi")
            g.load()
            # run_pipegraph()

        with RaisesDirectOrInsidePipegraph(ValueError):
            inner_tss()
        if ppg.util.global_pipegraph is not None:
            both_ppg_and_no_ppg.new_pipegraph()
        with RaisesDirectOrInsidePipegraph(ValueError):
            inner_tes()
        if ppg.util.global_pipegraph is not None:
            both_ppg_and_no_ppg.new_pipegraph()
        with RaisesDirectOrInsidePipegraph(ValueError):
            inner_chr()
        if ppg.util.global_pipegraph is not None:
            both_ppg_and_no_ppg.new_pipegraph()

    def test_alternative_loading_raises_on_missing_name(self):
        df = pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 5000,
                    "tes": 5500,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        )

        with pytest.raises(ValueError):
            genes.Genes(get_genome(), lambda: df)

    def test_alternative_loading_raises_on_invalid_chromosome(self):
        df = pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1b",
                    "strand": 1,
                    "tss": 5000,
                    "tes": 5500,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        )

        with RaisesDirectOrInsidePipegraph(ValueError):
            g = genes.Genes(get_genome(), lambda: df, name="shu")
            force_load(g.load())

    def test_alternative_loading_raises_on_non_int_tss(self):
        df = pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 5000.5,
                    "tes": 5500,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        )

        with RaisesDirectOrInsidePipegraph(ValueError):
            g = genes.Genes(get_genome(), lambda: df, name="shu")
            force_load(g.load())

    def test_alternative_loading_raises_on_non_int_tes(self):
        df = pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 5000,
                    "tes": "",
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        )

        with RaisesDirectOrInsidePipegraph(ValueError):
            g = genes.Genes(get_genome(), lambda: df, name="shu")
            force_load(g.load())

    def test_do_load_only_happens_once(self):
        df = pd.DataFrame(
            [
                {
                    "gene_stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 5000,
                    "tes": 5500,
                    "description": "bla",
                }
            ]
        )
        counter = [0]

        def load():
            counter[0] += 1
            return df

        g = genes.Genes(get_genome_chr_length(), load, name="shu")
        if ppg.inside_ppg():
            assert counter[0] == 0
            g.load()
            assert counter[0] == 0
            g.load()
            assert counter[0] == 0
            ppg.run_pipegraph()
        else:
            assert counter[0] == 1
            g.load()
            assert counter[0] == 1

    def test_filtering_away_works(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        g = genes.Genes(genome)
        filtered = g.filter("nogenes", lambda df: df["chr"] == "4")
        force_load(filtered.load())
        run_pipegraph()
        assert len(filtered.df) == 0
        assert "start" in filtered.df.columns
        assert "stop" in filtered.df.columns
        assert "tss" in filtered.df.columns
        assert "tes" in filtered.df.columns
        assert "gene_stable_id" in filtered.df.columns

    def test_annotators_are_kept_on_filtering(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        g = genes.Genes(genome)
        ca = Constant("shu", 5)
        g.add_annotator(ca)
        filtered = g.filter("nogenes", lambda df: df["chr"] == "4")
        assert filtered.has_annotator(ca)

    def test_filtering_returns_genes(self):
        g = genes.Genes(get_genome())
        on_chr_1 = g.filter("on_1", lambda df: df["chr"] == "1")
        assert g.__class__ == on_chr_1.__class__

    def test_overlap_genes_requires_two_genes(self):
        genome = get_genome()
        a = genes.Genes(genome)

        def sample_data():
            return pd.DataFrame(
                {"chr": ["Chromosome"], "start": [1000], "stop": [1100]}
            )

        b = regions.GenomicRegions("sha", sample_data, [], genome)
        force_load(a.load())
        force_load(b.load())
        run_pipegraph()

        with pytest.raises(ValueError):
            a.overlap_genes(b)

    def test_overlap_genes_raises_on_unequal_genomes(self):
        genome = get_genome("A")
        genomeB = get_genome("B")
        a = genes.Genes(genome)
        b = genes.Genes(genomeB)

        with pytest.raises(ValueError):
            a.overlap_genes(b)

    def test_overlap(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        g = genes.Genes(genome)
        on_chr_1 = g.filter("on_1", lambda df: df["chr"] == "1")
        on_chr_2 = g.filter("on_2", lambda df: df["chr"] == "2")
        one = g.filter("one", lambda df: df["gene_stable_id"] == "fake1")
        force_load(on_chr_1.load())
        force_load(on_chr_2.load())
        force_load(one.load())
        run_pipegraph()
        assert len(on_chr_1.df) == 2
        assert len(on_chr_2.df) == 1
        assert len(one.df) == 1
        assert g.overlap_genes(on_chr_1) == len(on_chr_1.df)
        assert on_chr_1.overlap_genes(g) == len(on_chr_1.df)
        assert on_chr_1.overlap_genes(on_chr_1) == len(on_chr_1.df)
        assert g.overlap_genes(on_chr_2) == len(on_chr_2.df)
        assert on_chr_2.overlap_genes(g) == len(on_chr_2.df)
        assert on_chr_2.overlap_genes(on_chr_2) == len(on_chr_2.df)
        assert g.overlap_genes(one) == len(one.df)
        assert one.overlap_genes(g) == len(one.df)
        assert one.overlap_genes(one) == len(one.df)

        assert on_chr_1.overlap_genes(one) == 1
        assert one.overlap_genes(on_chr_1) == 1

        assert on_chr_1.overlap_genes(on_chr_2) == 0
        assert on_chr_2.overlap_genes(on_chr_1) == 0

    def test_get_tss_regions(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        g = genes.Genes(genome)
        tss = g.regions_tss()
        force_load(tss.load())
        run_pipegraph()
        assert len(tss.df) == 3
        assert (tss.df["start"] == [5000, 5400, 5400]).all()
        assert (tss.df["stop"] == tss.df["start"] + 1).all()
        assert (tss.df["chr"] == ["1", "1", "2"]).all()

    def test_get_tes_regions(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 3000,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        g = genes.Genes(genome)
        tes = g.regions_tes()
        force_load(tes.load())
        run_pipegraph()
        assert len(tes.df) == 2
        assert (tes.df["start"] == [4900, 4900]).all()
        assert (tes.df["stop"] == tes.df["start"] + 1).all()
        assert (tes.df["chr"] == ["1", "2"]).all()

    def test_get_exons_regions_overlapping(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 3000,
                        "tes": 4900,
                        "description": "bla",
                        "name": "bla1",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                        "name": "bla2",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                        "name": "bla3",
                    },
                ]
            ),
            # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
            df_transcripts=pd.DataFrame(
                {
                    "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                    "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                    "chr": ["1", "1", "1", "2"],
                    "strand": [1, 1, -1, -1],
                    "start": [3100, 3000, 4910, 4900],
                    "stop": [4900, 4000, 5400, 5400],
                    "exons": [
                        [(3100, 4900)],
                        [(3000, 3500), (3300, 3330), (3750, 4000)],
                        [(4910, 5000), (5100, 5400)],
                        [(4900, 5400)],
                    ],
                }
            ),
        )
        g = genes.Genes(genome)
        exons = g.regions_exons_overlapping()
        force_load(exons.load())
        run_pipegraph()
        assert (exons.df["start"] == [3000, 3100, 3300, 3750, 4910, 5100, 4900]).all()
        assert (exons.df["stop"] == [3500, 4900, 3330, 4000, 5000, 5400, 5400]).all()
        assert (exons.df["chr"] == np.array(["1", "1", "1", "1", "1", "1", "2"])).all()

    def test_get_exons_regions_merging(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 3000,
                        "tes": 4900,
                        "description": "bla",
                        "name": "bla1",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                        "name": "bla2",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                        "name": "bla3",
                    },
                ]
            ),
            # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
            df_transcripts=pd.DataFrame(
                {
                    "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                    "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                    "chr": ["1", "1", "1", "2"],
                    "strand": [1, 1, -1, -1],
                    "start": [3100, 3000, 4910, 4900],
                    "stop": [4900, 4000, 5400, 5400],
                    "exons": [
                        [(3100, 4900)],
                        [(3000, 3500), (3300, 3330), (3750, 4000)],
                        [(4910, 5000), (5100, 5400)],
                        [(4900, 5400)],
                    ],
                }
            ),
        )
        g = genes.Genes(genome)
        exons = g.regions_exons_merged()
        force_load(exons.load())
        run_pipegraph()
        assert (exons.df["start"] == [3000, 4910, 5100, 4900]).all()
        assert (exons.df["stop"] == [4900, 5000, 5400, 5400]).all()
        assert (exons.df["chr"] == ["1", "1", "1", "2"]).all()

    def test_get_intron_regions(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 3000,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            ),
            # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
            df_transcripts=pd.DataFrame(
                {
                    "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                    "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                    "chr": ["1", "1", "1", "2"],
                    "strand": [1, 1, -1, -1],
                    "start": [3100, 3000, 4900, 4900],
                    "stop": [4900, 4000, 5400, 5400],
                    "exons": [
                        [(3100, 4900)],
                        [(3000, 3500), (3750, 4000)],
                        [(4900, 5000), (5100, 5400)],
                        [(4900, 5400)],
                    ],
                }
            ),
        )
        g = genes.Genes(genome)
        introns = g.regions_introns()
        force_load(introns.load())
        run_pipegraph()
        assert (introns.df["start"] == [3000, 3500, 4000, 5000]).all()
        assert (introns.df["stop"] == [3100, 3750, 4900, 5100]).all()
        # no intronic region on chr 2
        assert (introns.df["chr"] == ["1", "1", "1", "1"]).all()

     def test_get_gene_bodies(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 3000,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            ),
            # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
            df_transcripts=pd.DataFrame(
                {
                    "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                    "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                    "chr": ["1", "1", "1", "2"],
                    "strand": [1, 1, -1, -1],
                    "start": [3100, 3000, 4900, 4900],
                    "stop": [4900, 4000, 5400, 5400],
                    "exons": [
                        [(3100, 4900)],
                        [(3000, 3500), (3750, 4000)],
                        [(4900, 5000), (5100, 5400)],
                        [(4900, 5400)],
                    ],
                }
            ),
        )
        g = genes.Genes(genome)
        body = g.regions_body()
        force_load(body.load())
        run_pipegraph()
        assert (body.df["start"] == [3000, 4900, 4900]).all()
        assert (body.df["stop"] == [4901, 5401, 54001]).all()
        # no intronic region on chr 2
        assert (introns.df["chr"] == ["1", "1", "2"]).all()

@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenes:
    def test_write_bed(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        g = genes.Genes(genome)
        sample_filename = "genes.bed"
        g.write_bed(sample_filename)
        run_pipegraph()
        assert len(g.df) > 0
        read = read_bed(g.result_dir / sample_filename)
        assert len(read) == len(g.df)
        assert read[0].refseq == "1"
        assert read[1].refseq == "1"
        assert read[2].refseq == "2"
        assert read[0].position == 4900
        assert read[1].position == 5000
        assert read[2].position == 4900
        assert read[0].length == 500
        assert read[1].length == 500
        assert read[2].length == 500
        assert read[0].name == "fake2"
        assert read[1].name == "fake1"
        assert read[2].name == "fake3"

    def test_write_bed_auto_filename(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        g = genes.Genes(genome, sheet_name="da_genes")
        assert "/da_genes/" in str(g.result_dir)
        sample_filename = g.write_bed()[1]
        run_pipegraph()
        assert len(g.df) > 0
        read = read_bed(sample_filename)
        assert len(read) == len(g.df)
        assert read[0].refseq == "1"
        assert read[1].refseq == "1"
        assert read[2].refseq == "2"
        assert read[0].position == 4900
        assert read[1].position == 5000
        assert read[2].position == 4900
        assert read[0].length == 500
        assert read[1].length == 500
        assert read[2].length == 500
        assert read[0].name == "fake2"
        assert read[1].name == "fake1"
        assert read[2].name == "fake3"

    # def test_annotation_keeps_row_names(self):
    # g = genes.Genes(dummyGenome)
    # g.do_load()
    # row_names = g.df.row_names
    # g.annotate()
    # run_pipegraph()
    # self.assertTrue((row_names == g.df.row_names).all())

    def test_write(self):
        g = genes.Genes(get_genome())
        with pytest.raises(ValueError):
            g.write(mangler_function=lambda df: df.tail())
        a = g.write()
        b = g.write("b.xls")
        mangle = lambda df: df.head()  # noqa: E731
        c = g.write("c.xls", mangle)
        # this is ok... in ppg1 only., not in ppg2 strict
        # c = g.write("c.xls", mangle)
        if ppg.util.inside_ppg():  # this is ok outside of ppg
            if hasattr(ppg, "is_ppg2"):
                with pytest.raises(ValueError):
                    c = g.write("c.xls", mangle)
            with pytest.raises(ValueError):
                g.write("c.xls", lambda df: df.tail())
        run_pipegraph()
        afn = a[1]
        bfn = b[1]
        cfn = c[1]

        assert Path(afn).exists()
        assert Path(bfn).exists()
        assert Path(cfn).exists()
        assert_frame_equal(pd.read_csv(afn, sep="\t"), pd.read_excel(bfn))
        assert_frame_equal(
            pd.read_excel(bfn).head(),
            pd.read_excel(cfn),
            check_column_type=False,
            check_dtype=False,
        )

    def test_write_filtered(self):
        g = genes.Genes(get_genome())
        g2 = g.filter("filtered", lambda df: df.index[:2])
        g2.write(Path("filtered.xls").absolute())
        run_pipegraph()
        assert Path("filtered.xls").exists()
        df = pd.read_excel("filtered.xls")
        assert len(df) == 2
        assert "parent_row" in df.columns
        assert (df["parent_row"] == [0, 1]).all()

    def test_invalid_chromosomes(self):
        def a():
            return pd.DataFrame(
                {
                    "chr": "7a",
                    "start": 100,
                    "stop": 1000,
                    "tss": 100,
                    "tes": 1000,
                    "strand": 1,
                    "name": "gene1",
                    "gene_stable_id": "gene1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            genes.Genes(
                genome, alternative_load_func=a, name="my_genes", result_dir="my_genes"
            ).load()

    def test_invalid_tss(self):
        def a():
            return pd.DataFrame(
                {
                    "chr": "Chromosome",
                    "tss": "100",
                    "tes": 1000,
                    "strand": 1,
                    "name": "gene1",
                    "gene_stable_id": "gene1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            genes.Genes(
                genome, alternative_load_func=a, name="my_genes", result_dir="my_genes"
            ).load()

    def test_invalid_tes(self):
        def a():
            return pd.DataFrame(
                {
                    "chr": "Chromosome",
                    "tss": 100,
                    "tes": 1000.5,
                    "strand": 1,
                    "name": "gene1",
                    "gene_stable_id": "gene1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            genes.Genes(
                genome, alternative_load_func=a, name="my_genes", result_dir="my_genes"
            ).load()

    def test_invalid_start_stop(self):
        def a():
            return pd.DataFrame(
                {
                    "chr": "Chromosome",
                    "tss": 100,
                    "tes": 10,
                    "start": 100,
                    "stop": 10,
                    "strand": 1,
                    "name": "gene1",
                    "gene_stable_id": "gene1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            genes.Genes(
                genome, alternative_load_func=a, name="my_genes", result_dir="my_genes"
            ).load()


@pytest.mark.usefixtures("new_pipegraph")
class TestGenesPPG:
    def test_def_twice_alternative_loading_func(self):
        def a():
            return pd.DataFrame(
                {
                    "chr": "1",
                    "start": 100,
                    "stop": 1000,
                    "tss": 100,
                    "tes": 1000,
                    "strand": 1,
                    "name": "gene1",
                    "gene_stable_id": "gene1",
                },
                index=["gene1"],
            )

        def b():
            return pd.DataFrame(
                {
                    "chr": "1",
                    "start": 110,
                    "stop": 1000,
                    "tss": 110,
                    "tes": 1000,
                    "strand": 1,
                    "name": "gene1",
                    "gene_stable_id": "gene1",
                },
                index=["gene1"],
            )

        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    }
                ]
            )
        )
        gA = genes.Genes(
            genome, alternative_load_func=a, name="my_genes", result_dir="my_genes"
        )
        assert gA.result_dir.resolve() == Path("my_genes").resolve()
        gA.load()
        gA.load()
        with pytest.raises(ValueError):
            genes.Genes(genome, alternative_load_func=b, name="my_genes")


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenesFrom:
    def test_difference(self):
        genome = get_genome()
        a = genes.Genes(genome)
        b = a.filter("filtered", lambda df: df.index[:5])
        c = genes.FromDifference("delta", a, b)
        force_load(c.load())
        run_pipegraph()
        assert len(c.df) == len(a.df) - len(b.df)

    def test_intersection(self):
        genome = get_genome()
        a = genes.Genes(genome)
        b = a.filter("filtered", lambda df: df.index[:5])
        c = a.filter("filtered2", lambda df: df.index[4:6])
        with pytest.raises(ValueError):
            d = genes.FromIntersection("delta", b, c)
        d = genes.FromIntersection("delta", [b, c])
        force_load(a.load())
        force_load(d.load())
        run_pipegraph()
        assert len(d.df) == 1
        assert list(d.df.gene_stable_id) == list(a.df.gene_stable_id.loc[4:4])

    def test_intersection2(self):
        genome = get_genome()
        a = genes.Genes(genome)
        b = a.filter("filtered", lambda df: df.index[:5], vid="AA")
        c = b.filter("filtered2", lambda df: df.index[:1], vid=["BB", "CC"])
        with pytest.raises(ValueError):
            d = genes.FromIntersection("delta", b, c)
        d = genes.FromIntersection("delta", [b, c])
        force_load(a.load())
        force_load(d.load())
        run_pipegraph()
        assert len(d.df) == 1
        assert list(d.df.gene_stable_id) == list(a.df.gene_stable_id.loc[0:0])
        assert "AA" in d.vid
        assert "BB" in d.vid
        assert "CC" in d.vid

    def test_from_any(self):
        genome = get_genome()
        a = genes.Genes(genome)
        b = a.filter("filtered", lambda df: df.index[:5])
        c = a.filter("filtered2", lambda df: df.index[-5:])
        d = a.filter("filtered3", lambda df: df.index[10:15])
        e = genes.FromAny("delta", [b, c, d], sheet_name="shu")
        force_load(e.load())
        force_load(a.load())
        run_pipegraph()
        assert len(e.df) == 15
        assert sorted(list(e.df.gene_stable_id)) == sorted(
            list(a.df.gene_stable_id.iloc[:5])
            + list(a.df.gene_stable_id.iloc[10:15])
            + list(a.df.gene_stable_id.iloc[-5:])
        )
        assert "/shu/" in str(e.result_dir)

    def test_from_all(self):
        genome = get_genome()
        a = genes.Genes(genome)
        b = a.filter("filtered", lambda df: df.index[:5])
        c = a.filter("filtered2", lambda df: df.index[0:10])
        d = a.filter("filtered3", lambda df: df.index[3:10])
        e = genes.FromAll("delta", [b, c, d])
        force_load(e.load())
        force_load(a.load())
        run_pipegraph()
        assert len(e.df) == 2
        assert list(e.df.gene_stable_id) == list(a.df.gene_stable_id.loc[3:4])

    def test_from_none(self):
        genome = get_genome()
        a = genes.Genes(genome)
        b = a.filter("filtered", lambda df: df.index[:5])
        c = a.filter("filtered2", lambda df: df.index[-5:])
        d = a.filter("filtered3", lambda df: df.index[3:10])
        e = genes.FromNone("delta", [b, c, d])
        force_load(e.load())
        force_load(a.load())
        run_pipegraph()
        assert len(e.df) == len(a.df) - 5 - 5 - 5

    def test_genes_from_file(self, both_ppg_and_no_ppg):
        genome = get_genome()
        a = genes.Genes(genome)
        b = a.filter("filtered", lambda df: df.index[:5])
        b.write(Path("filtered.xls").absolute())
        force_load(b.load())
        print(both_ppg_and_no_ppg)
        run_pipegraph()
        assert not "summit middle" in a.df.columns
        assert not "summit middle" in b.df.columns
        print(both_ppg_and_no_ppg)
        both_ppg_and_no_ppg.new_pipegraph()
        genome = get_genome()
        c = genes.FromFile("reimport", genome, Path("filtered.xls").absolute())
        force_load(c.load())
        run_pipegraph()
        assert_frame_equal(b.df, c.df)

    def test_genes_from_file_of_transcripts(self):
        genome = get_genome()
        df = pd.DataFrame({"a column!": genome.df_transcripts.index[:5]})
        df.to_excel("transcripts.xlsx")
        a = genes.FromFileOfTranscripts(
            "my genes", genome, "transcripts.xlsx", "a column!"
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 5
        tr = set()
        for gene_stable_id in a.df["gene_stable_id"]:
            tr.update(
                [
                    tr.transcript_stable_id
                    for tr in genome.genes[gene_stable_id].transcripts
                ]
            )
        assert tr == set(genome.df_transcripts.index[:5])

    def test_genes_from_biotypes(self):
        genome = get_genome()
        nc = ["tRNA", "rRNA"]
        non_coding = genes.FromBiotypes(genome, nc)
        force_load(non_coding.load())
        run_pipegraph()
        assert len(non_coding.df) == genome.df_genes.biotype.isin(nc).sum()
