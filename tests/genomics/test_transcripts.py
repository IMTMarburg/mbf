import pypipegraph as ppg
import pandas as pd
import pytest
from pathlib import Path
import mbf_genomics.transcripts as transcripts
from mbf_genomics.annotator import Constant

from .shared import (
    get_genome,
    get_genome_chr_length,
    force_load,
    run_pipegraph,
    RaisesDirectOrInsidePipegraph,
    MockGenome,
)


@pytest.mark.usefixtures("new_pipegraph")
class TestTranscriptsLoadingPPGOnly:
    def test_loading_from_genome_is_singletonic(self):
        genome = get_genome()
        print(genome)
        genesA = transcripts.Transcripts(genome)
        genesB = transcripts.Transcripts(genome)
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
class TestTranscriptsLoading:
    def test_basic_loading_from_genome(self):
        g = transcripts.Transcripts(get_genome())
        force_load(g.load())
        run_pipegraph()
        assert len(g.df) == 246
        print(g.df)
        assert (g.df["gene_stable_id"][:3] == ["CRP_001", "CRP_002", "CRP_003"]).all()
        assert g.df["gene_stable_id"].iloc[-1] == "CRP_182"
        assert (
            g.df["transcript_stable_id"][:3] == ["BAF35032", "BAF35033", "BAF35034"]
        ).all()
        assert g.df["transcript_stable_id"].iloc[-1] == "BAF35213"
        assert g.df["start"].iloc[-1] == 158_648
        assert g.df["stop"].iloc[-1] == 159_662
        assert g.df["strand"].iloc[-1] == -1
        assert g.df["tss"].iloc[-1] == 159_662
        assert g.df["tes"].iloc[-1] == 158_648

    def test_filtering_with_annotator(self):
        import mbf_genomics

        mbf_genomics.transcripts.transcripts._transcripts_per_genome_singletons.clear()

        g = transcripts.Transcripts(get_genome())

        class CopyAnnoTTL(mbf_genomics.annotator.Annotator):
            def __init__(self):
                self.columns = ["copyAA"]

            def calc(self, df):
                return pd.DataFrame({"copyAA": df["transcript_stable_id"]})

        g += CopyAnnoTTL()
        filtered = g.filter("ax", ("transcript_stable_id", "==", "BAF35034"))
        force_load(filtered.annotate())
        run_pipegraph()
        print(filtered.name)
        print(g.name)
        print(filtered.df.columns)
        assert (filtered.df["gene_stable_id"] == ["CRP_003"]).all()
        print(filtered.df)
        assert (filtered.df["copyAA"] == ["BAF35034"]).all()

    def test_alternative_loading_raises_on_non_df(self):
        with RaisesDirectOrInsidePipegraph(ValueError):
            g = transcripts.Transcripts(get_genome_chr_length(), lambda: None, "myname")
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
            g = transcripts.Transcripts(get_genome(), lambda: df2, name="sha")
            g.load()
            # run_pipegraph()

        def inner_chr():
            df2 = df.copy()
            df2 = df2.drop("chr", axis=1)
            g = transcripts.Transcripts(get_genome(), lambda: df2, name="shu")
            g.load()
            # run_pipegraph()

        def inner_tes():
            df2 = df.copy()
            df2 = df2.drop("tes", axis=1)
            g = transcripts.Transcripts(get_genome(), lambda: df2, name="shi")
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
            transcripts.Transcripts(get_genome(), lambda: df)

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
            g = transcripts.Transcripts(get_genome(), lambda: df, name="shu")
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
            g = transcripts.Transcripts(get_genome(), lambda: df, name="shu")
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
            g = transcripts.Transcripts(get_genome(), lambda: df, name="shu")
            force_load(g.load())

    def test_do_load_only_happens_once(self):
        df = pd.DataFrame(
            [
                {
                    "gene_stable_id": "fake1",
                    "transcript_stable_id": "fake1",
                    "name": "fake1",
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

        g = transcripts.Transcripts(get_genome_chr_length(), load, name="shu")
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
            ),
            df_transcripts=pd.DataFrame(
                [
                    {
                        "stable_id": "fake1ts",
                        "gene_stable_id": "fake1",
                        "chr": "1",
                        "name": "fake1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                        "exons": [(5000, 5500)],
                    }
                ]
            ),
        )
        g = transcripts.Transcripts(genome)
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
            ),
            df_transcripts=pd.DataFrame(
                [
                    {
                        "stable_id": "fake1ts",
                        "gene_stable_id": "fake1",
                        "chr": "1",
                        "name": "fake1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                        "exons": [(5000, 5500)],
                    }
                ]
            ),
        )
        g = transcripts.Transcripts(genome)
        ca = Constant("shu", 5)
        g.add_annotator(ca)
        filtered = g.filter("nogenes", lambda df: df["chr"] == "4")
        assert filtered.has_annotator(ca)

    def test_filtering_returns_transcripts(self):
        g = transcripts.Transcripts(get_genome())
        on_chr_1 = g.filter("on_1", lambda df: df["chr"] == "1")
        assert g.__class__ == on_chr_1.__class__

    def test_get_tss_regions(self):
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
        g = transcripts.Transcripts(genome)
        tss = g.regions_tss()
        force_load(tss.load())
        run_pipegraph()
        assert len(tss.df) == 4
        assert (tss.df["start"] == [3000, 3100, 5400, 5400]).all()
        assert (tss.df["stop"] == tss.df["start"] + 1).all()
        assert (tss.df["chr"] == ["1", "1", "1", "2"]).all()
        assert (
            tss.df["transcript_stable_id"] == ["trans1b", "trans1a", "trans2", "trans3"]
        ).all()

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
        g = transcripts.Transcripts(genome)
        tes = g.regions_tes()
        force_load(tes.load())
        run_pipegraph()
        print(tes.df)
        assert len(tes.df) == 4
        assert (tes.df["start"] == [4000, 4900, 4900, 4900]).all()
        assert (tes.df["stop"] == tes.df["start"] + 1).all()
        assert (tes.df["chr"] == ["1", "1", "1", "2"]).all()
        assert (
            tes.df["transcript_stable_id"] == ["trans1b", "trans1a", "trans2", "trans3"]
        ).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenes:
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
                    "transcript_stable_id": "ts1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            transcripts.Transcripts(
                genome,
                alternative_load_func=a,
                name="my_genes",
                result_dir="my_genes",
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
                    "transcript_stable_id": "ts1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            transcripts.Transcripts(
                genome,
                alternative_load_func=a,
                name="my_genes",
                result_dir="my_genes",
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
                    "transcript_stable_id": "ts1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            transcripts.Transcripts(
                genome,
                alternative_load_func=a,
                name="my_genes",
                result_dir="my_genes",
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
                    "transcript_stable_id": "ts1",
                },
                index=["gene1"],
            )

        genome = get_genome()
        with RaisesDirectOrInsidePipegraph(ValueError):
            transcripts.Transcripts(
                genome,
                alternative_load_func=a,
                name="my_genes",
                result_dir="my_genes",
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
                    "transcript_stable_id": "ts1",
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
                    "transcript_stable_id": "ts1",
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
            ),
            df_transcripts=pd.DataFrame(
                {
                    "transcript_stable_id": ["trans1a"],
                    "gene_stable_id": ["fake1"],
                    "chr": ["1"],
                    "strand": [
                        1,
                    ],
                    "start": [5000],
                    "stop": [5500],
                    "exons": [
                        [(5000, 5500)],
                    ],
                }
            ),
        )
        gA = transcripts.Transcripts(
            genome, alternative_load_func=a, name="my_genes", result_dir="my_genes"
        )
        assert gA.result_dir.resolve() == Path("my_genes").resolve()
        gA.load()
        gA.load()
        with pytest.raises(ValueError):
            transcripts.Transcripts(genome, alternative_load_func=b, name="my_genes")
