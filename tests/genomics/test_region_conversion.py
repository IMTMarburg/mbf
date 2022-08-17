import pytest
import math
import pandas as pd
import pypipegraph as ppg

from mbf.genomics import regions, genes
from mbf.genomics.annotator import Constant

from .shared import (
    get_genome,
    get_genome_chr_length,
    force_load,
    run_pipegraph,
    RaisesDirectOrInsidePipegraph,
)


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenomicRegionConvertTests:
    def test_random_same_number(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1"],
                    "start": [10, 100, 1000],
                    "stop": [12, 110, 1110],
                    "column_that_will_disappear": ["A", "b", "c"],
                }
            )

        def convert(df):
            res = df[["chr", "start", "stop"]]
            res = res.assign(start=res["start"] + 1)
            return res

        if ppg.inside_ppg():
            deps = [ppg.ParameterInvariant("shuParam", ("hello"))]
        else:
            deps = []
        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        a.add_annotator(Constant("Constant", 5))
        a.annotate()
        b = a.convert("a+1", convert, dependencies=deps)
        force_load(b.load())
        for d in deps:
            assert d in b.load().lfg.prerequisites
        run_pipegraph()
        assert len(a.df) == len(b.df)
        assert (a.df["start"] == b.df["start"] - 1).all()
        assert "column_that_will_disappear" in a.df.columns
        assert not ("column_that_will_disappear" in b.df.columns)

    def test_raises_on_conversion_function_not_returning_df(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1"],
                    "start": [10, 100, 1000],
                    "stop": [12, 110, 1110],
                }
            )

        def convert(df):
            return None

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        with RaisesDirectOrInsidePipegraph(ValueError):
            a.convert("a+1", convert)
            force_load(a.load())

    def test_raises_on_non_genome(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1"],
                    "start": [10, 100, 1000],
                    "stop": [12, 110, 1110],
                }
            )

        def convert(df):
            res = df[:]
            res["start"] += 1
            return res

        genome = get_genome_chr_length()
        a = regions.GenomicRegions("sharum", sample_data, [], genome)
        a.convert("a+1", convert, genome)

        with pytest.raises(ValueError):
            a.convert("a+1b", convert, "hello")

    def test_grow(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "2"],
                    "start": [10, 100, 1000],
                    "stop": [12, 110, 1110],
                }
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("agrown", regions.convert.grow(12))
        force_load(b.load())
        run_pipegraph()
        assert (b.df["start"] == [0, 100 - 12, 1000 - 12]).all()
        assert (b.df["stop"] == [24, 110 + 12, 1110 + 12]).all()

    def test_grow2(self):
        def sample_data():
            return pd.DataFrame(
                [
                    {"chr": "1", "start": 7774885, "stop": 7791673},
                    {"chr": "1", "start": 8286026, "stop": 8298500},
                    {"chr": "1", "start": 8323232, "stop": 8342008},
                ]
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("agrown", regions.convert.grow(50000), on_overlap="merge")
        force_load(b.load())
        run_pipegraph()
        assert (b.df["start"] == [7724885, 8236026]).all()
        assert (b.df["stop"] == [7841673, 8392008]).all()

    def test_promotorize(self):
        import numpy as np

        g = genes.Genes(get_genome())
        b = g.convert("b", regions.convert.promotorize(444), on_overlap="ignore")
        c = g.convert("c", regions.convert.promotorize(444, 10), on_overlap="ignore")
        force_load(b.load())
        force_load(c.load())
        run_pipegraph()
        assert len(g.df) > 0
        assert len(g.df) == len(b.df) + 1  # we drop one that ends up at 0..0
        assert len(g.df) == len(c.df)  # guess it's now 0..10?
        assert "strand" in b.df.columns
        assert "strand" in c.df.columns
        # we have to go by index - the order might change
        # convert to list of strings - bug in at, it won't work otherwise
        b_df = b.df.assign(gene_stable_id=[x for x in b.df.gene_stable_id])
        g_df = g.df.assign(gene_stable_id=[x for x in g.df.gene_stable_id])
        c_df = c.df.assign(gene_stable_id=[x for x in c.df.gene_stable_id])
        b_df = b_df.set_index("gene_stable_id")
        c_df = c_df.set_index("gene_stable_id")
        g_df = g_df.set_index("gene_stable_id")
        assert set(b_df.index) == set(g_df[1:].index)  # again the one that we dropped
        assert set(c_df.index) == set(g_df.index)  #

        for ii in b_df.index:
            assert g_df.at[ii, "strand"] == b_df.at[ii, "strand"]
            assert g_df.at[ii, "strand"] == c_df.at[ii, "strand"]

            if g_df.at[ii, "strand"] == 1:
                assert b_df.at[ii, "start"] == max(0, g_df.at[ii, "tss"] - 444)
                assert b_df.at[ii, "stop"] == max(0, g_df.at[ii, "tss"])

                assert c_df.at[ii, "start"] == max(0, g_df.at[ii, "tss"] - 444)
                assert c_df.at[ii, "stop"] == max(0, g_df.at[ii, "tss"] + 10)
            else:
                assert b_df.at[ii, "start"] == max(0, g_df.at[ii, "tss"])
                assert b_df.at[ii, "stop"] == max(0, g_df.at[ii, "tss"] + 444)

                assert c_df.at[ii, "start"] == max(0, g_df.at[ii, "tss"] - 10)
                assert c_df.at[ii, "stop"] == max(0, g_df.at[ii, "tss"] + 444)
        assert c_df.start.dtype == np.int32
        assert b_df.start.dtype == np.int32
        assert c_df.stop.dtype == np.int32
        assert b_df.stop.dtype == np.int32

    def test_merge_connected(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": [
                        "1",
                        "1",
                        "1",
                        "1",
                        "1",
                        "2",
                        "3",
                        "3",
                        "3",
                        "4",
                        "4",
                        "4",
                        "5",
                    ],
                    "start": [
                        10,
                        13,
                        110,
                        300,
                        400,
                        102,
                        5,
                        6,
                        6000,
                        10,
                        100,
                        200,
                        100,
                    ],
                    "stop": [
                        18,
                        100,
                        200,
                        400,
                        410,
                        1000,
                        5000,
                        4900,
                        6010,
                        100,
                        150,
                        300,
                        110,
                    ],
                }
            )

        a = regions.GenomicRegions(
            "sharum", sample_data, [], get_genome_chr_length(), on_overlap="ignore"
        )
        b = a.convert("agrown", regions.convert.merge_connected())
        force_load(b.load())
        run_pipegraph()
        assert (b.df["chr"] == ["1", "1", "1", "2", "3", "3", "4", "4", "5"]).all()
        assert (b.df["start"] == [10, 110, 300, 102, 5, 6000, 10, 200, 100]).all()
        assert (b.df["stop"] == [100, 200, 410, 1000, 5000, 6010, 150, 300, 110]).all()

    def test_merge_connected_2(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "2"],
                    "start": [27897200, 27898600, 0],
                    "stop": [27897300, 27898700, 100],
                }
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("agrown", regions.convert.merge_connected())
        force_load(b.load())
        run_pipegraph()
        assert (b.df["start"] == [27897200, 27898600, 0]).all()
        assert (b.df["stop"] == [27897300, 27898700, 100]).all()

    def test_liftover(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "2"],
                    "start": [27897200, 27898600, 100000],
                    "stop": [27897300, 27898700, 100100],
                }
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("hg38", regions.convert.lift_over("hg19ToHg38"))
        force_load(b.load())
        run_pipegraph()
        # made these with the ucsc web liftover utility
        # http://genome.ucsc.edu/cgi-bin/hgLiftOver
        assert (b.df["start"] == [27570689, 27572089, 100000]).all()
        assert (b.df["stop"] == [27570789, 27572189, 100100]).all()

    def test_liftover_filter_chr(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["2", "1", "11_gl000202_random", "MT"],
                    "start": [100, 27897200, 500, 100000],
                    "stop": [10000, 27897300, 5000, 100100],
                    "copy": ["D", "A", "B", "C"],
                    "name": ["d", "a", "b", "c"],
                }
            )

        a = regions.GenomicRegions(
            "sharum",
            sample_data,
            [],
            get_genome_chr_length(
                {
                    "1": 100000,
                    "2": 100000,
                    "11_gl000202_random": 100000,
                    "MT": 100000,
                    "11": 1000000,
                }
            ),
        )
        b = a.convert(
            "hg38",
            regions.convert.lift_over(
                "hg19ToHg38", keep_name=True, filter_to_these_chromosomes=["1"]
            ),
        )
        force_load(b.load())
        run_pipegraph()
        # made these with the ucsc web liftover utility
        # http://genome.ucsc.edu/cgi-bin/hgLiftOver
        print(b.df)
        assert (b.df["start"] == [27570689]).all()
        assert (b.df["stop"] == [27570789]).all()
        assert (b.df["copy"] == ["A"]).all()
        assert (b.df["name"] == ["a"]).all()
        assert (b.df["chr"] == ["1"]).all()

    def test_windows(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["1", "2"], "start": [0, 0], "stop": [200, 99 * 3]}
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("ndiwo", regions.convert.windows(99, False))
        c = a.convert("ndiwo2", regions.convert.windows(99, True))
        force_load(b.load())
        force_load(c.load())
        run_pipegraph()
        assert (b.df["start"] == [0, 99, 99 * 2, 0, 99, 99 * 2]).all()
        assert (b.df["stop"] == [99, 99 * 2, 200, 99, 99 * 2, 99 * 3]).all()
        assert (c.df["start"] == [0, 99, 0, 99, 99 * 2]).all()
        assert (c.df["stop"] == [99, 99 * 2, 99, 99 * 2, 99 * 3]).all()

    def test_cookiecutter_summit(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["1", "2"], "start": [0, 0], "stop": [200, 99 * 3]}
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("cookie", regions.convert.cookie_summit(a.summit_annotator, 220))
        c = a.convert(
            "cookieB",
            regions.convert.cookie_summit(
                a.summit_annotator, 220, drop_those_outside_chromosomes=True
            ),
        )
        force_load(b.load())
        force_load(c.load())
        run_pipegraph()
        assert len(b.df) == 2
        assert len(c.df) == 1
        assert (b.df["start"] == [0, math.floor(99 * 3 / 2) - 110]).all()
        assert (b.df["stop"] == [100 + 110, math.floor(99 * 3 / 2) + 110]).all()

    def test_name_must_be_string(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["1", "2"], "start": [0, 0], "stop": [200, 99 * 3]}
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        with pytest.raises(ValueError):
            a.convert(123, regions.convert.shift(50))

    def test_shift(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["1", "2"], "start": [0, 0], "stop": [200, 99 * 3]}
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("cookie", regions.convert.shift(50))
        force_load(b.load())
        run_pipegraph()
        assert len(b.df) == 2
        assert (b.df["start"] == [50, 50]).all()
        assert (b.df["stop"] == [200 + 50, 99 * 3 + 50]).all()

    def test_summit(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["1", "2"], "start": [0, 0], "stop": [200, 99 * 3]}
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("cookie", regions.convert.summit(a.summit_annotator))
        force_load(b.load())
        run_pipegraph()
        assert len(b.df) == 2
        assert (b.df["start"] == [100, 148]).all()
        assert (b.df["stop"] == [101, 149]).all()

    def test_cookie_cutter(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2"],
                    "start": [0, 0],
                    "stop": [200, 99 * 3],
                    "strand": [-1, 1],
                }
            )

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome_chr_length())
        b = a.convert("cookie", regions.convert.cookie_cutter(100))
        force_load(b.load())
        run_pipegraph()
        assert (b.df["start"] == [50, 98]).all()
        assert (b.df["stop"] == [50 + 100, 98 + 100]).all()
        assert (b.df["strand"] == [-1, 1]).all()
