import logging

logging.basicConfig(level=logging.WARNING)
import pytest
import numpy
import os
import pypipegraph as ppg
import pandas as pd
from pathlib import Path
from pandas.testing import assert_frame_equal
import dppd
import dppd_plotnine  # noqa:F401
from mbf.qualitycontrol.testing import assert_image_equal
from mbf.sampledata import get_sample_data


import mbf.genomics.regions as regions
from mbf.genomics.annotator import Constant, Annotator

from .shared import (
    get_genome,
    get_genome_chr_length,
    force_load,
    inside_ppg,
    run_pipegraph,
    RaisesDirectOrInsidePipegraph,
    MockGenome,
)

dp, X = dppd.dppd()


@pytest.mark.usefixtures("new_pipegraph")
class TestGenomicRegionsLoadingPPGOnly:
    def test_dependency_passing(self):
        job = ppg.ParameterInvariant("sha", (None,))
        a = regions.GenomicRegions("shu", lambda: None, [job], get_genome())
        load_job = a.load()
        assert job in load_job.lfg.prerequisites

    def test_dependency_may_be_iterable_instead_of_list(self):
        job = ppg.ParameterInvariant("shu", (None,))
        a = regions.GenomicRegions("shu", lambda: None, (job,), get_genome())
        load_job = a.load()
        assert job in load_job.lfg.prerequisites

    def test_depenencies_must_be_jobs(self):
        ppg.ParameterInvariant("shu", (None,))
        if hasattr(ppg, "is_ppg2"):
            expected = KeyError
        else:
            expected = ValueError
        with pytest.raises(expected):  # ppg2 KeyError
            regions.GenomicRegions("shu", lambda: None, ["shu"], get_genome())


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenomicRegionsLoading:
    def test_raises_on_duplicate_name(self, both_ppg_and_no_ppg):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["Chromosome"], "start": [1000], "stop": [1100]}
            )

        regions.GenomicRegions("shu", sample_data, [], get_genome())

        if inside_ppg():
            with pytest.raises(ValueError):
                regions.GenomicRegions("shu", sample_data, [], get_genome())
            both_ppg_and_no_ppg.new_pipegraph()
            regions.GenomicRegions(
                "shu", sample_data, [], get_genome()
            )  # should not raise

    def test_raises_on_non_iterable_dependencies(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["Chromosome"], "start": [1000], "stop": [1100]}
            )

        with pytest.raises(ValueError):
            regions.GenomicRegions("shu", sample_data, "aaeu", get_genome())
        with pytest.raises(ValueError):
            regions.GenomicRegions("shu", sample_data, 1, get_genome())
        with pytest.raises(ValueError):
            regions.GenomicRegions("shu", sample_data, iter([]), get_genome())

    def test_loading(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["Chromosome"], "start": [1000], "stop": [1100]}
            )

        a = regions.GenomicRegions("sha", sample_data, [], get_genome())
        if inside_ppg():
            assert not hasattr(a, "df")
            force_load(a.load())
        else:
            assert hasattr(a, "df")
        run_pipegraph()
        assert hasattr(a, "df")
        assert len(a.df) == 1
        assert "chr" in a.df.columns
        assert "start" in a.df.columns
        assert "stop" in a.df.columns

    def test_filtering_copy_anno(self, clear_annotators):
        import mbf.genomics

        def sample_data():
            return pd.DataFrame(
                {
                    "chr": "Chromosome",
                    "start": [1000, 1001, 1002],
                    "stop": [1100, 1101, 1102],
                }
            )

        a = regions.GenomicRegions(
            "sha", sample_data, [], get_genome(), on_overlap="ignore"
        )
        b = a.filter("filtered", ("start", "==", 1001))

        class CopyAnno(mbf.genomics.annotator.Annotator):
            def __init__(self):
                self.columns = ["copy"]

            def calc(self, df):
                return pd.DataFrame({"copy": df["start"]})

        a += CopyAnno()

        if inside_ppg():
            assert not hasattr(a, "df")
            force_load(a.load())
            force_load(b.annotate())
        else:
            assert hasattr(a, "df")
        run_pipegraph()
        print(b.df)
        assert (b.df["start"] == [1001]).all()
        assert (b.df["copy"] == [1001]).all()

    def test_raises_on_invalid_on_overlap(self):
        def inner():
            regions.GenomicRegions(
                "shu",
                lambda: None,
                [],
                get_genome(),
                on_overlap="run in circles all about",
            )

        with pytest.raises(ValueError):
            inner()

    def test_magic(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["Chromosome"], "start": [1000], "stop": [1100]}
            )

        a = regions.GenomicRegions("shu", sample_data, [], get_genome())
        hash(a)
        str(a)
        repr(a)
        bool(a)
        a.load()
        run_pipegraph()
        with pytest.raises(TypeError):
            iter(a)

    def test_loading_missing_start(self):
        def sample_data():
            return pd.DataFrame({"chr": "1", "stop": [1100]})

        g = get_genome()

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sha", sample_data, [], g)
            force_load(a.load)

    def test_loading_missing_chr(self):
        def sample_data():
            return pd.DataFrame({"start": [1000], "stop": [1100]})

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sha", sample_data, [], get_genome())
            force_load(a.load)

    def test_loading_missing_stop(self):
        def sample_data():
            return pd.DataFrame({"chr": "Chromosome", "start": [1200]})

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sha", sample_data, [], get_genome())
            force_load(a.load)

    def test_loading_raises_on_invalid_chromosome(self):
        def sample_data():
            return pd.DataFrame({"chr": ["1b"], "start": [1200], "stop": [1232]})

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
            force_load(a.load)

    def test_loading_raises_on_no_int_start(self):
        def sample_data():
            return pd.DataFrame(
                {"chr": ["Chromosome"], "start": ["shu"], "stop": [1232]}
            )

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
            force_load(a.load)

    def test_loading_raises_on_no_int_stop(self):
        def sample_data():
            return pd.DataFrame({"chr": ["Chromosome"], "start": [2], "stop": [20.0]})

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
            force_load(a.load)

    def test_loading_raises_on_no_str_chr(self):
        def sample_data():
            return pd.DataFrame({"chr": [1], "start": [2], "stop": [20]})

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
            force_load(a.load)

    def test_loading_raises_on_not_dataframe(self):
        def sample_data():
            return None

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
            force_load(a.load)

    def test_loading_raises_on_overlapping(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["Chromosome", "Chromosome"],
                    "start": [1000, 1010],
                    "stop": [1100, 1020],
                }
            )

        with RaisesDirectOrInsidePipegraph(ValueError):

            a = regions.GenomicRegions(
                "sha", sample_data, [], get_genome(), on_overlap="raise"
            )
            force_load(a.load)

    def test_raises_on_negative_intervals(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": [
                        "Chromosome",
                        "Chromosome",
                        "Chromosome",
                        "Chromosome",
                        "Chromosome",
                    ],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [8, 110, 1110, 11110, 111_110],
                }
            )

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
            force_load(a.load)

    def test_raises_on_overlapping_intervals(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                }
            )

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
            force_load(a.load)

    def test_index_reset(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["Chromosome"],
                    "start": [1000],
                    "stop": [1100],
                    "myindex": ["a"],
                }
            ).set_index("myindex")

        a = regions.GenomicRegions("sharum", sample_data, [], get_genome())
        force_load(a.load)
        run_pipegraph()
        assert a.df.index == ["a"]
        assert not "myindex" in a.df.columns

    def test_merges_overlapping_intervals(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="merge"
        )
        force_load(a.load)
        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 1110

    def test_merge_overlapping_with_function(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                    "pick_me": [1, 2, 3, 4, 5],
                }
            )

        def merge_function(subset_df):
            row = subset_df.iloc[0].to_dict()
            row["pick_me"] = numpy.max(subset_df["pick_me"])
            return row

        a = regions.GenomicRegions(
            "shu",
            sample_data,
            [],
            get_genome_chr_length(),
            on_overlap=("merge", merge_function),
        )
        force_load(a.load)
        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 1110
        assert a.df.iloc[0]["pick_me"] == 3

    def test_merge_overlapping_with_function2(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                    "pick_me": [10, 2, 3, 4, 5],
                }
            )

        def merge_function(subset_df):
            row = subset_df.iloc[0].to_dict()
            row["pick_me"] = numpy.max(subset_df["pick_me"])
            return row

        a = regions.GenomicRegions(
            "shu",
            sample_data,
            [],
            get_genome_chr_length(),
            on_overlap=("merge", merge_function),
        )
        force_load(a.load)
        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 1110
        assert a.df.iloc[0]["pick_me"] == 10

    def test_merge_overlapping_with_function_ignores_returned_intervals(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                    "pick_me": [1, 2, 3, 4, 5],
                }
            )

        def merge_function(subset_df):
            row = subset_df.iloc[0].to_dict()
            row["start"] = 9000
            row["pick_me"] = numpy.max(subset_df["pick_me"])
            return row

        a = regions.GenomicRegions(
            "shu",
            sample_data,
            [],
            get_genome_chr_length(),
            on_overlap=("merge", merge_function),
        )
        force_load(a.load)
        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 1110
        assert a.df.iloc[0]["pick_me"] == 3

    def test_merge_overlapping_with_function_raises_on_non_dict(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                    "pick_me": [1, 2, 3, 4, 5],
                }
            )

        def merge_function(subset_df):
            row = subset_df.iloc[0].copy()
            row["pick_me"] = numpy.max(subset_df["pick_me"])
            return row

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions(
                "shu",
                sample_data,
                [],
                get_genome_chr_length(),
                on_overlap=("merge", merge_function),
            )
            force_load(a.load)

    def test_merge_overlapping_with_function_raises_on_unknown_column(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                    "pick_me": [1, 2, 3, 4, 5],
                }
            )

        def merge_function(subset_df):
            row = subset_df.iloc[0][:]
            row["does not exist"] = row["pick_me"]
            return row

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions(
                "shu",
                sample_data,
                [],
                get_genome_chr_length(),
                on_overlap=("merge", merge_function),
            )
            force_load(a.load)

    def test_merge_overlapping_with_function_raises_on_missing(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                    "pick_me": [1, 2, 3, 4, 5],
                }
            )

        def merge_function(subset_df):
            row = subset_df.iloc[0]
            del row["pick_me"]
            return None

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions(
                "shu",
                sample_data,
                [],
                get_genome_chr_length(),
                on_overlap=("merge", merge_function),
            )
            force_load(a.load)

    def test_merges_overlapping_intervals_next_to_each_other(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["Chromosome", "Chromosome"],
                    "start": [10, 21],
                    "stop": [20, 100],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome(), on_overlap="merge"
        )
        force_load(a.load)
        run_pipegraph()
        assert (a.df["start"] == [10, 21]).all()
        assert (a.df["stop"] == [20, 100]).all()

    def test_merges_overlapping_earlier_overlaps_later_with_one_in_between(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["Chromosome", "Chromosome", "Chromosome", "Chromosome"],
                    "start": [3100, 3000, 3750, 4910],
                    "stop": [4900, 3500, 4000, 5000],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome(), on_overlap="merge"
        )
        force_load(a.load)

        run_pipegraph()
        assert (a.df["start"] == [3000, 4910]).all()
        assert (a.df["stop"] == [4900, 5000]).all()

    def test_merges_overlapping_intervals_multiple(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5", "1"],
                    "start": [10, 100, 1000, 10000, 100_000, 1005],
                    "stop": [1010, 110, 1110, 11110, 111_110, 2000],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="merge"
        )
        force_load(a.load)

        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 2000

    def test_merges_overlapping_intervals_multiple_with_function(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5", "1"],
                    "start": [10, 100, 1000, 10000, 100_000, 1005],
                    "stop": [1010, 110, 1110, 11110, 111_110, 2000],
                    "pick_me": [23, 1234, 2, 4, 5, 50],
                }
            )

        def merge_function(subset_df):
            row = subset_df.iloc[0].to_dict()
            row["pick_me"] = numpy.max(subset_df["pick_me"])
            return row

        a = regions.GenomicRegions(
            "shu",
            sample_data,
            [],
            get_genome_chr_length(),
            on_overlap=("merge", merge_function),
        )
        force_load(a.load)

        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 2000
        assert a.df.iloc[0]["pick_me"] == 50

    def test_on_overlap_drop(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5", "1", "5", "5", "5"],
                    "start": [10, 100, 1000, 10000, 100_000, 1005, 5000, 6049, 6000],
                    "stop": [1010, 110, 1110, 11110, 111_110, 2000, 5050, 6051, 6050],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="drop"
        )
        force_load(a.load)
        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 100
        assert a.df.iloc[0]["stop"] == 110
        assert a.df.iloc[1]["start"] == 10000
        assert a.df.iloc[1]["stop"] == 11110
        assert a.df.iloc[2]["start"] == 5000
        assert a.df.iloc[2]["stop"] == 5050
        assert a.df.iloc[3]["start"] == 100_000
        assert a.df.iloc[3]["stop"] == 111_110

    def test_on_overlap_drop_nested(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5", "1", "5", "5", "5"],
                    "start": [10, 100, 1000, 10000, 100_000, 1005, 5000, 6049, 6000],
                    "stop": [1010, 110, 1110, 11110, 111_110, 2000, 5050, 6100, 6050],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="drop"
        )
        force_load(a.load)
        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 100
        assert a.df.iloc[0]["stop"] == 110
        assert a.df.iloc[1]["start"] == 10000
        assert a.df.iloc[1]["stop"] == 11110
        assert a.df.iloc[2]["start"] == 5000
        assert a.df.iloc[2]["stop"] == 5050
        assert a.df.iloc[3]["start"] == 100_000
        assert a.df.iloc[3]["stop"] == 111_110

    def test_on_overlap_drop_empty(self):
        def sample_data():
            return pd.DataFrame({"chr": [], "start": [], "stop": []})

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="drop"
        )
        force_load(a.load)
        run_pipegraph()
        assert len(a.df) == 0

    def test_on_overlap_ignore(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5", "1", "5"],
                    "start": [10, 100, 1000, 10000, 100_000, 1005, 5000],
                    "stop": [1010, 110, 1110, 11110, 111_110, 2000, 5050],
                    "index": ["a", "b", "c", "d", "e", "f", "g"],
                    "is_overlapping": 55,
                }
            ).set_index("index")

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="ignore"
        )
        force_load(a.load)

        run_pipegraph()
        assert len(a.df) == 7
        assert (a.df.index == ["a", "c", "f", "b", "d", "g", "e"]).all()
        assert not (a.df.is_overlapping == 55).any()  # make sure we drop the column
        assert (
            a.df.is_overlapping == [True, True, False, False, False, False, False]
        ).any()

    def test_merging_apperantly_creating_negative_intervals(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "start": [
                        140_688_139,
                        140_733_871,
                        140_773_241,
                        140_792_022,
                        141_032_547,
                        141_045_565,
                        141_069_938,
                        141_075_938,
                        141_098_775,
                        141_108_518,
                        141_131_159,
                        -4423,
                        -4352,
                        -3398,
                        -3329,
                        -1770,
                        -1693,
                        -737,
                        3400,
                        -598,
                    ],
                    "chr": [
                        str(x)
                        for x in [
                            "1",
                            "1",
                            "1",
                            "1",
                            "1",
                            "1",
                            "1",
                            "1",
                            "1",
                            "1",
                            "1",
                            "2",
                            "2",
                            "2",
                            "2",
                            "2",
                            "2",
                            "2",
                            "2",
                            "2",
                        ]
                    ],
                    "stop": [
                        140_767_241,
                        140_786_022,
                        141_026_547,
                        141_039_565,
                        141_064_437,
                        141_070_437,
                        141_092_775,
                        141_102_518,
                        141_125_159,
                        141_145_045,
                        141_213_431,
                        1577,
                        1648,
                        2602,
                        2671,
                        4230,
                        4307,
                        5263,
                        9400,
                        5402,
                    ],
                }
            )

        with RaisesDirectOrInsidePipegraph(
            ValueError, "All starts need to be positive"
        ):
            a = regions.GenomicRegions(
                "shu", sample_data, [], get_genome_chr_length(), on_overlap="merge"
            )
            force_load(a.load)

    def test_merge_test(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "1", "3", "5"],
                    "start": [10, 10, 1000, 10000, 100_000],
                    "stop": [100, 100, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="merge"
        )
        force_load(a.load)

        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 100

    def test_merge_identical_ok(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "1", "3", "5"],
                    "start": [10, 10, 1000, 10, 100_000],
                    "stop": [100, 100, 1110, 100, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu",
            sample_data,
            [],
            get_genome_chr_length(),
            on_overlap="merge_identical",
        )
        force_load(a.load)

        run_pipegraph()
        assert len(a.df) == 4
        assert a.df.iloc[0]["start"] == 10
        assert a.df.iloc[0]["stop"] == 100

    def test_merge_identical_raises(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "1", "3", "5"],
                    "start": [10, 10, 1000, 10000, 100_000],
                    "stop": [100, 120, 1110, 11110, 111_110],
                }
            )

        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions(
                "shu",
                sample_data,
                [],
                get_genome_chr_length(),
                on_overlap="merge_identical",
            )
            force_load(a.load)

    def test_regions_merge_in_init_does_not_add_strand_if_it_was_missing(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "1", "3", "5"],
                    "start": [10, 10, 1000, 10000, 100_000],
                    "stop": [100, 120, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="merge"
        )
        force_load(a.load())
        run_pipegraph()
        assert not "strand" in a.df.columns

    def test_plot_job(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                }
            )

        x = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="merge"
        )
        pj = x.plot(
            "shu.png",
            lambda df: dp(df)
            .p9()
            .add_scatter("chr", "start")
            .pd,  # also tests the write-to-result_dir_part
        )
        fn = "results/GenomicRegions/shu/shu.png"
        if inside_ppg():
            assert isinstance(pj[0], ppg.FileGeneratingJob)
            assert pj[1].absolute() == Path(fn).absolute()
        else:
            assert str(pj[1]) == str(Path(fn).absolute())
        run_pipegraph()
        assert_image_equal(fn)

    def test_plot_job_with_custom_calc_function(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                }
            )

        x = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(), on_overlap="merge"
        )

        def calc(df):
            df = df.copy()
            df = df.assign(shu=["shu"] * len(df))
            return df

        pj = x.plot(
            Path("shu.png").absolute(),
            lambda df: dp(df).p9().add_scatter("shu", "start"),
            calc,
        )
        fn = str(Path("shu.png").absolute())
        if inside_ppg():
            assert isinstance(pj[0], ppg.FileGeneratingJob)
            assert str(pj[1].absolute()) == fn
        else:
            assert str(pj[1]) == fn
        run_pipegraph()
        assert_image_equal(fn)


@pytest.mark.usefixtures("new_pipegraph")
class TestGenomicRegionsAnnotationDependencyies:
    def setUp(self):
        def sample_data():
            df = pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [12, 110, 1110, 11110, 111_110],
                }
            )
            df = df.assign(summit=df["start"] + (df["stop"] - df["start"]) / 2)
            return df

        self.a = regions.GenomicRegions("shu", sample_data, [], get_genome_chr_length())

    def test_anno_job_depends_on_load(self):
        self.setUp()
        if inside_ppg():
            assert not hasattr(self.a, "df")
        else:
            assert hasattr(self.a, "df")
        ca = Constant("Constant", 5)
        anno_job = self.a.add_annotator(ca)
        assert isinstance(anno_job(), ppg.Job)


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenomicRegionsAnnotation:
    def setUp(self):
        def sample_data():
            df = pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [12, 110, 1110, 11110, 111_110],
                }
            )
            df = df.assign(summit=df["start"] + (df["stop"] - df["start"]) / 2)
            return df

        self.a = regions.GenomicRegions("shu", sample_data, [], get_genome_chr_length())

    def test_anno_jobs_are_singletonic(self):
        self.setUp()
        ca = Constant("Constant", 5)
        assert len(self.a.annotators) == 1
        anno_job = self.a.add_annotator(ca)
        anno_job2 = self.a.add_annotator(ca)
        assert anno_job is anno_job2

    def test_anno_jobs_are_singletonic_across_names(self):
        self.setUp()
        ca = Constant("Constant", 5)
        assert len(self.a.annotators) == 1
        anno_job = self.a.add_annotator(ca)
        cb = Constant("Constant", 5)
        assert ca is cb
        anno_job2 = self.a.add_annotator(cb)
        assert anno_job is anno_job2

    def test_has_annotator(self):
        self.setUp()
        ca = Constant("Constant", 5)
        assert not self.a.has_annotator(ca)
        self.a.add_annotator(ca)
        assert self.a.has_annotator(ca)

    def test_annotator_by_name(self):
        self.setUp()
        ca = Constant("Constant", 5)
        assert not self.a.has_annotator(ca)
        self.a.add_annotator(ca)
        assert ca == self.a.get_annotator(ca.columns[0])

    def test_anno_jobs_add_columns(self):
        self.setUp()
        ca = Constant("Constant", 5)
        assert len(self.a.annotators) == 1
        self.a.add_annotator(ca)
        force_load(self.a.annotate(), "test_anno_jobs_add_columns")
        run_pipegraph()
        assert ca.columns[0] in self.a.df.columns


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestGenomicRegionsWriting:
    def setUp(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                    "name": ["a", "b", "c", "d", "e"],
                    "notname": ["A", "B", "C", "D", "E"],
                    "strand": [1, 1, 1, 1, -1],
                }
            )

        self.a = regions.GenomicRegions("shu", sample_data, [], get_genome_chr_length())
        self.sample_filename = str(Path("sample.dat").absolute())
        try:
            os.unlink(self.sample_filename)
        except OSError:
            pass

    def test_write_bed(self):
        self.setUp()
        from mbf.fileformats.bed import read_bed

        self.a.write_bed(self.sample_filename)
        run_pipegraph()
        assert len(self.a.df) > 0
        read = read_bed(self.sample_filename)
        assert len(read) == len(self.a.df)
        assert read[0].refseq == "1"
        assert read[1].refseq == "1"
        assert read[2].refseq == "2"
        assert read[3].refseq == "3"
        assert read[4].refseq == "5"
        assert read[0].position == 10
        assert read[1].position == 1000
        assert read[2].position == 100
        assert read[3].position == 10000
        assert read[4].position == 100_000
        assert read[0].length == 1
        assert read[1].length == 110
        assert read[2].length == 10
        assert read[3].length == 1110
        assert read[4].length == 11110
        assert read[0].name == "Noname"

    def test_write_bed_with_name(self):
        self.setUp()
        from mbf.fileformats.bed import read_bed

        self.a.write_bed(self.sample_filename, region_name="name")
        run_pipegraph()
        assert len(self.a.df) > 0
        read = read_bed(self.sample_filename)
        assert len(read) == len(self.a.df)
        assert read[0].refseq == "1"
        assert read[1].refseq == "1"
        assert read[2].refseq == "2"
        assert read[3].refseq == "3"
        assert read[4].refseq == "5"
        assert read[0].position == 10
        assert read[1].position == 1000
        assert read[2].position == 100
        assert read[3].position == 10000
        assert read[4].position == 100_000
        assert read[0].length == 1
        assert read[1].length == 110
        assert read[2].length == 10
        assert read[3].length == 1110
        assert read[4].length == 11110
        assert read[0].name == "a"
        assert read[1].name == "c"
        assert read[2].name == "b"
        assert read[3].name == "d"
        assert read[4].name == "e"

    def test_write_bigbed_name_column(self):
        self.setUp()
        from mbf.fileformats.bed import read_bigbed

        self.a.write_bigbed(self.sample_filename, "notname")
        run_pipegraph()
        assert len(self.a.df) > 0
        read = read_bigbed(
            self.sample_filename, {"1": 10000, "2": 20000, "3": 30000, "5": 500_000}
        )
        print(read)
        should = self.a.df.reset_index(drop=True)

        assert len(read) == len(should)
        assert (read["chr"] == should["chr"]).all()
        assert (read["start"] == should["start"]).all()
        assert (read["stop"] == should["stop"]).all()
        assert (read["strand"] == should["strand"]).all()
        assert (
            read["name"].str.upper() == should["name"].str.upper()
        ).all()  # bigbed seems to store uprcse names?

    def test_write_bigbed(self):
        self.setUp()
        from mbf.fileformats.bed import read_bigbed

        self.a.write_bigbed(self.sample_filename)
        run_pipegraph()
        assert len(self.a.df) > 0
        read = read_bigbed(
            self.sample_filename, {"1": 10000, "2": 20000, "3": 30000, "5": 500_000}
        )
        should = self.a.df.reset_index(drop=True)
        print(read)
        print(should)

        assert len(read) == len(should)
        assert (read["chr"] == should["chr"]).all()
        assert (read["start"] == should["start"]).all()
        assert (read["stop"] == should["stop"]).all()
        assert (read["strand"] == should["strand"]).all()

    def test_write_bed_with_name_column_not_found(self):
        self.setUp()

        with RaisesDirectOrInsidePipegraph(KeyError):
            self.a.write_bed(self.sample_filename, region_name="name_not_found")

    def test_write(self):
        self.setUp()
        self.a.write(self.sample_filename)
        run_pipegraph()
        assert len(self.a.df) > 0
        df = pd.read_csv(self.sample_filename, sep="\t")
        df["chr"] = df["chr"].astype(str)
        assert_frame_equal(
            df, self.a.df.reset_index(drop=True)
        )  # , check_less_precise=2)

    def test_write_without_filename(self):
        self.setUp()
        self.a.result_dir = Path("")
        self.a.write()
        run_pipegraph()
        assert os.path.exists("shu.tsv")
        os.unlink("shu.tsv")

    def test_write_sorted(self):
        self.setUp()
        # sorting by chromosome means they would've been equal anyhow, since we
        # internally sort by chr
        self.a.write(self.sample_filename, lambda df: df.sort_values("start"))
        run_pipegraph()
        assert len(self.a.df) > 0
        df = pd.read_csv(self.sample_filename, sep="\t")
        df["chr"] = df["chr"].astype(str)
        assert set(df.columns) == set(self.a.df.columns)
        df = df[self.a.df.columns]  # make them have the same order
        assert not (df["start"].values == self.a.df["start"].values).all()
        assert (df == self.a.df.sort_values("start").reset_index(drop=True)).all().all()

    def test_plot_plots(self):
        self.setUp()

        pj = self.a.plot("shu.png", lambda df: dp(df).p9().add_scatter("start", "stop"))
        fn = "results/GenomicRegions/shu/shu.png"
        if inside_ppg():
            assert isinstance(pj[0], ppg.FileGeneratingJob)
            assert Path(pj[0].filenames[0]).absolute() == Path(fn).absolute()
            assert Path(pj[1]).absolute() == Path(fn).absolute()
        else:
            assert str(pj[1]) == str(Path(fn).absolute())
        run_pipegraph()
        assert_image_equal(fn)


class GRNameAnnotator(Annotator):
    def __init__(self, column_name="GR"):
        self.columns = [column_name]

    def calc_ddf(self, genomic_regions):
        """Return a dataframe with new columns, len(new_df) == len(genomic_regions.df).
        Merging is done by the genomic_region"""
        return pd.DataFrame(
            {self.columns[0]: [genomic_regions.name] * len(genomic_regions.df)}
        )


@pytest.mark.usefixtures("new_pipegraph")
class TestFilterTestDependencies:
    def setUp(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        self.genome = get_genome_chr_length(genome=get_genome())
        self.a = regions.GenomicRegions("shu", sample_data, [], self.genome)

    def test_dependecies(self):
        self.setUp()
        job = ppg.ParameterInvariant("shu_param", (None,))
        b = self.a.filter("sha", lambda df: df["chr"] == "1", dependencies=[job])
        force_load(b.load())
        assert job in b.load().lfg.prerequisites
        run_pipegraph()
        assert len(b.df) == 2

    def test_filter_can_depend_on_anno_jobs(self):
        self.setUp()
        anno = Constant("Constant", 5)
        anno_job = self.a.add_annotator(anno)

        def filter(df):
            return df[anno.columns[0]] == 5

        b = self.a.filter("sha", filter, dependencies=[anno_job])
        b.write("shu.tsv")
        run_pipegraph()
        # self.assertTrue(len(b.df) == 5)

    def test_filter_can_depend_on_single_job(self):
        self.setUp()
        anno = Constant("Constant", 5)
        anno_job = self.a.add_annotator(anno)

        def filter(df):
            return df[anno.columns[0]] == 5

        b = self.a.filter("sha", filter, dependencies=[anno_job])
        b.write("shu.tsv")
        run_pipegraph()
        assert len(b.df) == 5


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestFilter:
    def setUp(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        self.genome = get_genome_chr_length(genome=get_genome())
        self.a = regions.GenomicRegions("shu", sample_data, [], self.genome)

    def test_filtering(self):
        self.setUp()
        b = self.a.filter("sha", lambda df: df["chr"] == "1")
        force_load(b.load)
        run_pipegraph()
        assert len(b.df) == 2
        assert (b.df["start"] == [10, 1000]).all()

    def test_filter_raises_on_no_params(self):
        self.setUp()
        with pytest.raises(TypeError):
            self.a.filter("no_parms")

    def test_sorting(self):
        self.setUp()
        b = self.a.filter("sha", lambda df: df.sort_values("chr").index)
        force_load(b.load)
        run_pipegraph()
        assert len(b.df) == len(self.a.df)
        assert (b.df["chr"] == ["1", "1", "2", "3", "5"]).all()
        assert (b.df["start"] == [10, 1000, 100, 10000, 100_000]).all()

    def test_select_top_k(self):
        self.setUp()
        b = self.a.filter(
            "sha", lambda df: df.sort_values("start", ascending=False)[:2].index
        )
        force_load(b.load)
        run_pipegraph()
        assert len(b.df) == 2
        assert b.df.iloc[1]["chr"] == "5"
        assert b.df.iloc[0]["chr"] == "3"
        assert b.df.iloc[1]["start"] == 100_000
        assert b.df.iloc[0]["start"] == 10000
        assert b.df.iloc[1]["stop"] == 111_110
        assert b.df.iloc[0]["stop"] == 11110

    def test_select_top_k_raises_on_non_int(self):
        self.setUp()
        with RaisesDirectOrInsidePipegraph(TypeError):
            b = self.a.filter(
                "sha", lambda df: df.sort_values("start", ascending=False)[:2.0]
            )
            force_load(b.load)

    def test_annotator_inheritance(self):
        self.setUp()
        anno = GRNameAnnotator()
        self.a.add_annotator(anno)
        b = self.a.filter("sha", lambda df: df["chr"] == "1")
        assert b.has_annotator(anno)
        force_load(self.a.annotate, "test_annotator_inheritance")
        force_load(b.annotate, "test_annotator_inheritanceb")
        run_pipegraph()
        assert len(self.a.df) == 5
        assert len(b.df) == 2
        assert self.a.df.iloc[0][anno.columns[0]] == self.a.name
        assert b.df.iloc[0][anno.columns[0]] == self.a.name

    def test_filter_remove_overlapping(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame({"chr": ["1"], "start": [1001], "stop": [1002]})

        b = regions.GenomicRegions("b", sample_data, [], self.genome)
        c = regions.GenomicRegions_FilterRemoveOverlapping("c", self.a, b)
        c.write("shu.tsv")
        run_pipegraph()
        assert len(c.df) == 4
        assert (c.df["start"] == [10, 100, 10000, 100_000]).all()

    def test_filter_to_overlapping(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame({"chr": ["1"], "start": [1001], "stop": [1002]})

        b = regions.GenomicRegions("b", sample_data, [], self.genome)
        c = regions.GenomicRegions_FilterToOverlapping("c", self.a, b)
        c.write("shu_c.tsv")
        d = regions.GenomicRegions_FilterToOverlapping("d", self.a, [b, b])
        d.write("shu_d.tsv")
        run_pipegraph()
        print("a", self.a.df)
        print("b", b.df)
        print("c", c.df)
        assert len(c.df) == 1
        assert (c.df["chr"] == ["1"]).all()
        assert (c.df["start"] == [1000]).all()
        assert (c.df["stop"] == [1110]).all()
        assert_frame_equal(c.df, d.df)

    def test_filtering_keeps_non_annotated_non_canonical_columns(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                    "hello": [1, 2, 3, 4, 5],
                }
            )

        b = regions.GenomicRegions("b", sample_data, [], get_genome_chr_length())
        c = b.filter("c", df_filter_function=lambda df: df["chr"] == "1")
        c.write("shu.tsv")
        run_pipegraph()
        assert len(c.df) == 2
        assert (c.df["chr"] == ["1", "1"]).all()
        assert (c.df["start"] == [10, 1000]).all()
        assert (c.df["hello"] == [1, 3]).all()


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestInterval:
    def setUp(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        self.genome = get_genome_chr_length(genome=get_genome())
        self.a = regions.GenomicRegions("shu", sample_data, [], self.genome)

    def test_overlapping(self):
        self.setUp()
        force_load(self.a.load())
        run_pipegraph()
        found = self.a.get_overlapping("1", 800, 1001)
        assert len(found) == 1

    def test_overlapping_excludes_right_border(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame({"chr": ["1"], "start": [100], "stop": [200]})

        a = regions.GenomicRegions("a", sample_data, [], self.genome)
        force_load(a.load())
        run_pipegraph()
        print(dir(a))
        assert not len(a.get_overlapping("1", 5, 99))
        assert not len(a.get_overlapping("1", 5, 100))
        assert len(a.get_overlapping("1", 5, 101))

    def test_end_exclusion(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame(
                [
                    {"chr": "1", "start": 7_724_885, "stop": 7_841_673},
                    {"chr": "1", "start": 8_236_026, "stop": 8_392_008},
                ]
            )

        a = regions.GenomicRegions("sha", sample_data, [], self.genome)
        force_load(a.load())
        run_pipegraph()
        # emtpy query
        assert not a.has_overlapping("1", 7_778_885, 7_778_885)
        assert a.has_overlapping("1", 7_778_885, 7_778_886)
        ov = a.get_overlapping("1", 7_778_885, 7_778_885)
        assert len(ov) == 0
        ov = a.get_overlapping("1", 7_778_885, 7_778_886)
        assert len(ov) == 1
        assert ov.iloc[0]["start"] == 7_724_885

    def test_has_overlapping(self):
        self.setUp()
        force_load(self.a.load())
        run_pipegraph()
        assert self.a.has_overlapping("1", 800, 1001)

    def test_overlapping_point(self):
        self.setUp()
        force_load(self.a.load())
        run_pipegraph()
        found = self.a.get_overlapping("2", 105, 106)
        assert len(found) == 1
        found = self.a.get_overlapping("2", 105, 105)
        assert len(found) == 0

    def test_has_overlapping_point(self):
        self.setUp()
        force_load(self.a.load())
        run_pipegraph()
        assert not self.a.has_overlapping("2", 105, 105)
        assert self.a.has_overlapping("2", 105, 106)

    def test_non_overlapping(self):
        self.setUp()
        force_load(self.a.load())
        run_pipegraph()
        found = self.a.get_overlapping("1", 1200, 1300)
        assert len(found) == 0
        assert isinstance(found, pd.DataFrame)

        found = self.a.get_overlapping("4", 1200, 1300)
        assert len(found) == 0
        assert isinstance(found, pd.DataFrame)

    def test_has_overlapping_not(self):
        self.setUp()
        force_load(self.a.load())
        run_pipegraph()
        assert not self.a.has_overlapping("1", 1200, 1300)

    def test_load_works_with_empty_df(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame({"chr": [], "start": [], "stop": []})

        b = regions.GenomicRegions("shub", sample_data, [], self.genome)
        force_load(b.load())
        run_pipegraph()
        assert len(b.df) == 0
        assert not len(b.get_overlapping("1", 0, 10000))
        assert not b.has_overlapping("1", 0, 10000)
        assert not len(b.get_closest_by_start("1", 500))

    def test_load_works_with_missing_chrs(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame({"chr": ["2"], "start": [100], "stop": [1000]})

        b = regions.GenomicRegions(
            "shub", sample_data, [], get_genome_chr_length(genome=self.genome)
        )
        force_load(b.load())
        run_pipegraph()
        assert len(b.df) == 1
        assert not len(b.get_overlapping("1", 0, 10000))
        assert not b.has_overlapping("1", 0, 10000)
        assert not len(b.get_closest_by_start("1", 500))

    def test_load_works_with_empty_df_handle_overlapping(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame({"chr": [], "start": [], "stop": []})

        b = regions.GenomicRegions(
            "shub", sample_data, [], self.genome, on_overlap="ignore"
        )
        force_load(b.load())
        run_pipegraph()
        assert len(b.df) == 0
        assert not len(b.get_overlapping("1", 0, 10000))
        assert not b.has_overlapping("1", 0, 10000)
        assert not len(b.get_closest_by_start("1", 500))

    def test_closest(self):
        self.setUp()
        force_load(self.a.load())
        run_pipegraph()
        found = self.a.get_closest_by_start("1", 5)
        assert len(found) == 1
        assert found.iloc[0]["chr"] == "1"
        assert found.iloc[0]["start"] == 10

        found = self.a.get_closest_by_start("2", 150)  # before, but not after...
        assert len(found) == 1
        assert found.iloc[0]["chr"] == "2"
        assert found.iloc[0]["start"] == 100

        # that's still closer to the left one..
        found = self.a.get_closest_by_start("1", 501)
        assert len(found) == 1
        assert found.iloc[0]["chr"] == "1"
        assert found.iloc[0]["start"] == 10

        # definatly closer to the second one...
        found = self.a.get_closest_by_start("1", 701)
        assert len(found) == 1
        assert found.iloc[0]["chr"] == "1"
        assert found.iloc[0]["start"] == 1000

        found = self.a.get_closest_by_start("1", 1050)  # within an interval...
        assert len(found) == 1
        assert found.iloc[0]["chr"] == "1"
        assert found.iloc[0]["start"] == 1000
        found = self.a.get_closest_by_start("1", 505)  # test the border
        assert len(found) == 1
        assert found.iloc[0]["chr"] == "1"
        assert found.iloc[0]["start"] == 10
        found = self.a.get_closest_by_start("1", 506)  # test the border
        assert len(found) == 1
        assert found.iloc[0]["chr"] == "1"
        assert found.iloc[0]["start"] == 1000

        found = self.a.get_closest_by_start("4", 701)  # empty chromosome...
        assert len(found) == 0


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestIntervalTestsNeedingOverlapHandling(TestInterval):
    def test_nested(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "1", "1", "3"],
                    "start": [10, 25, 30, 100, 0],
                    "stop": [25, 35, 35, 110, 100],
                }
            )

        genome = get_genome_chr_length(genome=get_genome())
        a = regions.GenomicRegions("shu", sample_data, [], genome, on_overlap="ignore")
        force_load(a.load())
        run_pipegraph()
        assert (a.get_overlapping("1", 25, 26).start == [25]).all()
        assert (a.get_overlapping("1", 15, 26).start == [10, 25]).all()
        assert (a.get_overlapping("1", 30, 40).start == [25, 30]).all()

    def test_closest_overlapping(self):
        self.setUp()

        def sample_data():
            return pd.DataFrame(
                {
                    "start": [
                        566_564,
                        569_592,
                        713_866,
                        935_162,
                        1_051_311,
                        1_279_151,
                        1_282_803,
                        1_310_387,
                        1_337_193,
                        1_447_089,
                    ],
                    "chr": [
                        "1",
                        "1",
                        "1",
                        "1",
                        "1",
                        "1",
                        "1",
                        "1",
                        "1",
                        "1",
                    ],
                    "stop": [
                        667_063,
                        570_304,
                        714_288,
                        937_142,
                        1_052_403,
                        1_281_233,
                        1_283_631,
                        1_311_060,
                        1_337_881,
                        1_447_626,
                    ],
                }
            )

        b = regions.GenomicRegions(
            "shbu", sample_data, [], self.genome, on_overlap="ignore"
        )
        force_load(b.load())
        run_pipegraph()
        closest = b.get_closest_by_start("1", 1_015_327)
        assert len(closest)
        assert closest.iloc[0]["chr"] == "1"
        assert closest.iloc[0]["start"] == 1_051_311
        closest = b.get_closest_by_start("1", 0)
        assert closest.iloc[0]["start"] == 566_564
        closest = b.get_closest_by_start("1", 566_564)
        assert closest.iloc[0]["start"] == 566_564
        closest = b.get_closest_by_start("1", 569_590)
        assert closest.iloc[0]["start"] == 569_592
        closest = b.get_closest_by_start("1", 569_592)
        assert len(closest) == 1
        assert closest.iloc[0]["start"] == 569_592
        closest = b.get_closest_by_start("1", 569_592 + 5)
        assert len(closest) == 1
        assert closest.iloc[0]["start"] == 569_592
        assert len(b.get_overlapping("1", 569_592, 569_593)) == 2
        closest = b.get_closest_by_start("1", 667_063)
        assert closest.iloc[0]["start"] == 713_866
        closest = b.get_closest_by_start("1", 570_304)
        assert closest.iloc[0]["start"] == 569_592

    def test_empty_chromosome(self):
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
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )

        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "1", "1", "1", "1"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        b = regions.GenomicRegions("shbu", sample_data, [], genome, on_overlap="ignore")
        force_load(b.load())
        run_pipegraph()
        closest = b.get_closest_by_start("2", 10)
        assert len(closest) == 0


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestAssortedGenomicRegion:
    def test_get_no_of_entries(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(genome=get_genome())
        )
        force_load(a.load)
        run_pipegraph()
        assert a.get_no_of_entries() == 5

    def test_get_covered_bases(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(genome=get_genome())
        )
        force_load(a.load)
        run_pipegraph()
        assert a.covered_bases == 1 + 10 + 110 + 1110 + 11110

    def test_get_covered_bases_raises_on_possible_overlaps(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu",
            sample_data,
            [],
            get_genome_chr_length(genome=get_genome()),
            on_overlap="ignore",
        )
        force_load(a.load())
        run_pipegraph()

        with pytest.raises(ValueError):
            a.covered_bases

    def test_get_mean_size(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(genome=get_genome())
        )
        force_load(a.load)
        run_pipegraph()
        assert a.mean_size == (1.0 + 10 + 110 + 1110 + 11110) / 5


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestSetOperationsOnGenomicRegions:
    def sample_to_gr(self, a, name, on_overlap="merge", genome=None):
        if not hasattr(self, "genome"):
            self.genome = get_genome()

        def sample_a(a=a):
            data = {"chr": [], "start": [], "stop": []}
            if genome is not None:
                chr = list(genome.get_chromosome_lengths())[0]
            else:
                chr = "Chromosome"
            for start, stop in a:
                data["chr"].append(chr)
                data["start"].append(start)
                data["stop"].append(stop)
            return pd.DataFrame(data)

        a = regions.GenomicRegions(
            name,
            sample_a,
            [],
            self.genome if genome is None else genome,
            on_overlap=on_overlap,
        )
        return a

    def handle(self, a, b, should, operation, multi_pass=False):
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        if multi_pass:
            c = operation("c", [a, b])
        else:
            c = operation("c", a, b)
        c.write("shu.tsv")
        run_pipegraph()
        print(c.df)
        assert len(should) == len(c.df)
        starts = []
        stops = []
        for start, stop in should:
            starts.append(start)
            stops.append(stop)
        if len(starts) == 0:
            assert len(c.df) == 0
        else:
            assert (c.df["start"].values == starts).all()
            assert (c.df["stop"].values == stops).all()
        return a, b, c

    def test_union(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        should = [(10, 120), (400, 450), (600, 700)]
        self.handle(a, b, should, regions.GenomicRegions_Union, multi_pass=True)

    def test_union_preserves_annos(self):
        a = [(10, 100), (400, 450)]
        ca1 = Constant("one", 1)
        b = [(80, 120), (600, 700)]
        ca2 = Constant("two", 2)
        a = self.sample_to_gr(a, "a")
        a.add_annotator(ca1)
        b = self.sample_to_gr(b, "b")
        b.add_annotator(ca2)
        c = regions.GenomicRegions_Union("c", [a, b])
        assert not c.has_annotator(ca1)
        assert not c.has_annotator(ca2)

    def test_difference(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        should = [(10, 80), (400, 450)]
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_difference_split_joins(self):
        a = [(0, 1000)]
        b = [(500, 1000)]
        should = [(0, 500)]
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_difference_removes_annos(self):
        a = [(10, 100), (400, 450)]
        ca1 = Constant("one", 1)
        b = [(80, 120), (600, 700)]
        ca2 = Constant("two", 2)
        a = self.sample_to_gr(a, "a")
        a.add_annotator(ca1)
        b = self.sample_to_gr(b, "b")
        b.add_annotator(ca2)
        c = regions.GenomicRegions_Difference("c", a, b)
        assert not c.has_annotator(ca1)
        assert not c.has_annotator(ca2)
        c.copy_annotators(a, b)
        assert c.has_annotator(ca1)
        assert c.has_annotator(ca2)

    def test_difference_start(self):
        a = [(10, 100), (400, 450)]
        b = [(10, 20), (600, 700)]
        should = [(20, 100), (400, 450)]
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_difference_start_2(self):
        a = [(10, 100), (400, 450)]
        b = [(5, 20), (600, 700)]
        should = [(20, 100), (400, 450)]
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_difference_split(self):
        a = [(100, 1000)]
        b = [(80, 120), (500, 600), (800, 1200)]
        should = [(120, 500), (600, 800)]
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_difference_split_2(self):
        a = [(10, 100), (400, 450)]
        b = [(11, 20), (600, 700)]
        should = [(10, 11), (20, 100), (400, 450)]
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_difference_empty(self):
        a = [(10, 100), (400, 450)]
        b = [(5, 2000)]
        should = []
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_difference_adjoining_regions(self):
        # in the end, this is a fence-post error.
        a = [(495, 687)]
        b = [(124, 495)]
        should = [(495, 687)]
        self.handle(a, b, should, regions.GenomicRegions_Difference)

    def test_intersection(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        should = [(80, 100)]
        self.handle(a, b, should, regions.GenomicRegions_Intersection)

    def test_from_common(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        c = regions.GenomicRegions_Common("c", [a, b])
        c.write()
        run_pipegraph()
        assert (c.df["start"] == [10]).all()
        assert (c.df["stop"] == [120]).all()

    def test_from_common_in_at_least(self):
        a = [(10, 100), (400, 450), (1000, 1200)]
        b = [(80, 120), (600, 700), (1000, 1100)]
        c = [(70, 110)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        c = self.sample_to_gr(c, "c")
        d = regions.GenomicRegions_CommonInAtLeastX("d", [a, b, c], 3)
        d.write()
        e = regions.GenomicRegions_CommonInAtLeastX("e", [a, b, c], 2)
        e.write()
        run_pipegraph()
        assert (d.df["start"] == [10]).all()
        assert (d.df["stop"] == [120]).all()
        assert (e.df["start"] == [10, 1000]).all()
        assert (e.df["stop"] == [120, 1200]).all()

    def test_intersection_preserves_annos(self):
        a = [(10, 100), (400, 450)]
        ca1 = Constant("one", 1)
        b = [(80, 120), (600, 700)]
        ca2 = Constant("two", 2)
        a = self.sample_to_gr(a, "a")
        a.add_annotator(ca1)
        b = self.sample_to_gr(b, "b")
        b.add_annotator(ca2)
        c = regions.GenomicRegions_Intersection("c", a, b)
        assert not c.has_annotator(ca1)
        assert not c.has_annotator(ca2)
        c.copy_annotators(a, b)
        assert c.has_annotator(ca1)
        assert c.has_annotator(ca2)

    def test_intersection_empty(self):
        a = [(10, 100), (400, 450)]
        b = [(600, 700)]
        should = []
        self.handle(a, b, should, regions.GenomicRegions_Intersection)

    def test_overlapping(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        should = [(10, 120)]
        self.handle(a, b, should, regions.GenomicRegions_Overlapping)

    def test_overlapping_empty(self):
        a = [(10, 100), (400, 450)]
        b = [(110, 120), (600, 700)]
        should = []
        self.handle(a, b, should, regions.GenomicRegions_Overlapping)

    def test_overlapping_preserves_annos(self):
        a = [(10, 100), (400, 450)]
        ca1 = Constant("one", 1)
        b = [(80, 120), (600, 700)]
        ca2 = Constant("two", 2)
        a = self.sample_to_gr(a, "a")
        a.add_annotator(ca1)
        b = self.sample_to_gr(b, "b")
        b.add_annotator(ca2)
        c = regions.GenomicRegions_Overlapping("c", a, b)
        assert not c.has_annotator(ca1)
        assert not c.has_annotator(ca2)
        c.copy_annotators(a, b)
        assert c.has_annotator(ca1)
        assert c.has_annotator(ca2)

    def test_invert(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1"],
                    "start": [10, 100, 1000],
                    "stop": [12, 110, 1110],
                }
            )

        g = get_genome_chr_length(genome=get_genome())
        a = regions.GenomicRegions("sharum", sample_data, [], g)
        b = regions.GenomicRegions_Invert("b", a)

        def should_data():
            return pd.DataFrame(
                {
                    "chr": [
                        "1",
                        "1",
                        "1",
                        "2",
                        "2",
                        "3",
                        "4",
                        "5",
                    ],  # interval-less chromosomes become totally covered...
                    "start": [0, 12, 1110, 0, 110, 0, 0, 0],
                    "stop": [
                        10,
                        1000,
                        100_000,
                        100,
                        200_000,
                        300_000,
                        400_000,
                        500_000,
                    ],
                }
            )

        should = regions.GenomicRegions("should", should_data, [], g)
        force_load(b.load())
        force_load(should.load())
        run_pipegraph()
        assert_frame_equal(
            b.df.reset_index(drop=True), should.df.reset_index(drop=True)
        )

    def test_invert_twice(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1"],
                    "start": [10, 100, 1000],
                    "stop": [12, 110, 1110],
                }
            )

        a = regions.GenomicRegions(
            "sharum", sample_data, [], get_genome_chr_length(genome=get_genome())
        )
        b = regions.GenomicRegions_Invert("b", a)
        c = regions.GenomicRegions_Invert("c", b)
        force_load(c.load())
        run_pipegraph()
        assert not "strand" in a.df.columns
        assert not "strand" in b.df.columns
        assert not "strand" in c.df.columns
        assert_frame_equal(
            c.df.reset_index(drop=True), a.df.reset_index(drop=True), check_dtype=False
        )

    def test_invert_preserves_annos(self):
        a = [(10, 100), (400, 450)]
        ca1 = Constant("one", 1)
        a = self.sample_to_gr(a, "a")
        a.add_annotator(ca1)
        c = regions.GenomicRegions_Invert("c", a)
        assert not c.has_annotator(ca1)
        c.copy_annotators(a)
        assert c.has_annotator(ca1)

    def test_overlap_basepairs(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        assert a.overlap_basepair(b) == 20
        assert b.overlap_basepair(a) == 20

    def test_overlap_basepairs_raise_on_overlaps(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        a = self.sample_to_gr(a, "a", on_overlap="ignore")
        b = self.sample_to_gr(b, "b")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        with pytest.raises(ValueError):
            a.overlap_basepair(b)

    def test_overlap_basepairs_identical(self):
        a = [(10, 100), (400, 450)]
        b = a
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        assert a.overlap_basepair(b) == a.covered_bases

    def test_overlap_basepairs_identical_plus_overlapping(self):
        a = [(10, 100), (400, 450), (498, 499)]
        b = [(10, 100), (400, 425), (412, 499)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b", on_overlap="ignore")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        assert a.overlap_basepair(b) == a.covered_bases

    def test_overlap_percentage(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        assert a.overlap_percentage(b) == 20 / (40 + 100.0)
        assert b.overlap_percentage(a) == 20 / (40 + 100.0)

    def test_overlap_percentage_one_empty(self):
        a = [(10, 100), (400, 450)]
        b = []
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        assert a.overlap_percentage(b) == 0
        assert b.overlap_percentage(a) == 0

    def test_intersection_count(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        assert a.intersection_count(b) == 1
        assert b.intersection_count(a) == 1

    def test_intersection_count_unequal_filter_count(self):

        a = [(500, 600)]
        b = [(400, 510), (590, 650)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        just_a = regions.GenomicRegions_FilterRemoveOverlapping("a-b", a, b)
        just_b = regions.GenomicRegions_FilterRemoveOverlapping("b-a", b, a)
        force_load(a.load())
        force_load(b.load())
        force_load(just_a.load())
        force_load(just_b.load())
        run_pipegraph()
        # coming this way, we only have one intersection
        assert len(a.df) - len(just_a.df) == a.intersection_count(b)
        # coming this way, it's two of them
        assert len(b.df) - len(just_b.df) != b.intersection_count(a)

    def test_overlap_count(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        assert a.overlap_count(b) == 1
        assert b.overlap_count(a) == 1

    def test_overlap_count_raises_on_unequal_genomes(self):
        a = [(10, 100), (400, 450)]
        b = [(80, 120), (600, 700)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(
            b,
            "b",
            genome=MockGenome(
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
            ),
        )
        a.load()
        force_load(a.load())
        b.load()
        force_load(b.load())
        run_pipegraph()
        with pytest.raises(ValueError):
            a.overlap_count(b)

    def test_overlap_count_slightly_more_complex(self):

        a = [(500, 600), (1000, 1100), (3000, 3001)]
        b = [(400, 510), (590, 650), (1050, 1055), (2000, 2002)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        force_load(a.load())
        force_load(b.load())
        just_a = regions.GenomicRegions_FilterRemoveOverlapping("just_a", a, b)
        just_b = regions.GenomicRegions_FilterRemoveOverlapping("just_b", b, [a])
        force_load(just_a.load())
        force_load(just_b.load())
        both = regions.GenomicRegions_Overlapping("ab", a, b)
        force_load(both.load())
        run_pipegraph()
        assert a.overlap_count(b) == 2
        assert b.overlap_count(a) == 2
        assert len(just_a.df) == 1
        assert len(just_b.df) == 1
        assert a.overlap_count(b) == len(both.df)
        assert b.overlap_count(a) == len(both.df)
        # self.assertEqual(len(a) + len(b), len(just_a) + len(just_b) + 2 * a.overlap_count(b))
        # self.assertEqual(len(a) + len(b), len(just_a) + len(just_b) + 2 * a.overlap_count(a))

    def test_overlap_count_extending(self):

        a = [(500, 600), (650, 1000)]
        b = [(400, 510), (590, 700), (10000, 10001)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        b.load()
        m1 = regions.GenomicRegions_Overlapping("ab", a, b)
        m2 = regions.GenomicRegions_Overlapping("ba", b, a)
        force_load(m1.load())
        force_load(m2.load())
        run_pipegraph()
        assert len(m1.df) == 1
        assert len(m2.df) == 1
        assert a.overlap_count(b) == 1
        assert b.overlap_count(a) == 1

    def test_overlap_count_extending_does_not_miss_last(self):

        a = [(500, 600), (650, 1000)]
        b = [(400, 510), (590, 700)]
        a = self.sample_to_gr(a, "a")
        b = self.sample_to_gr(b, "b")
        a.load()
        b.load()
        m1 = regions.GenomicRegions_Overlapping("ab", a, b)
        m2 = regions.GenomicRegions_Overlapping("ba", b, a)

        force_load(m1.load())
        force_load(m2.load())
        run_pipegraph()
        assert len(m1.df) == 1
        assert len(m2.df) == 1
        assert a.overlap_count(b) == 1
        assert b.overlap_count(a) == 1

    def test_raises_on_unequal_genome(self):
        def sample_data():
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [11, 110, 1110, 11110, 111_110],
                }
            )

        g = get_genome()
        a = regions.GenomicRegions(
            "shu", sample_data, [], get_genome_chr_length(genome=g)
        )
        genome2 = get_genome_chr_length(genome=get_genome(name="arosebyanyothername"))
        b = regions.GenomicRegions("sha", sample_data, [], genome2)

        def inner_union():
            regions.GenomicRegions_Union("shu", [a, b])

        def inner_intersection():
            regions.GenomicRegions_Intersection("sha", a, b)

        def inner_difference():
            regions.GenomicRegions_Difference("sha", a, b)

        def inner_overlapping():
            regions.GenomicRegions_Overlapping("sha", a, b)

        def inner_overlap_basepair():
            a.overlap_basepair(b)

        def inner_overlap_percentag():
            a.overlap_percentage(b)

        def inner_intersection_count():
            a.intersection_count(b)

        with pytest.raises(ValueError):
            inner_union()
        with pytest.raises(ValueError):
            inner_intersection()
        with pytest.raises(ValueError):
            inner_difference()
        with pytest.raises(ValueError):
            inner_overlapping()
        with pytest.raises(ValueError):
            inner_overlap_basepair()
        with pytest.raises(ValueError):
            inner_overlap_percentag()
        with pytest.raises(ValueError):
            inner_intersection_count()
        force_load(a.load())
        force_load(b.load())
        run_pipegraph()

        def inner_iter_intersections():
            next(a._iter_intersections(b))

        with pytest.raises(ValueError):
            inner_iter_intersections()


@pytest.mark.usefixtures("both_ppg_and_no_ppg")
class TestFromXYZ:
    def test_gff(self):
        a = regions.GenomicRegions_FromGFF(
            "shu",
            get_sample_data("mbf_genomics/test.gff3"),
            get_genome_chr_length(genome=get_genome()),
            filter_function=lambda entry: entry["source"]
            == b"Regions_of_sig_enrichment",
            comment_char="#",
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 4
        assert (a.df["chr"] == ["2", "3", "3", "3"]).all()
        assert (a.df["start"] == [662_743, 53968, 58681, 68187]).all()
        assert (
            a.df["score"]
            == [
                b"-0.583140831299777",
                b"0.41667200444801",
                b"1.20483507668127",
                b"0.742024408736128",
            ]
        ).all()

    def test_gff_with_name(self):
        a = regions.GenomicRegions_FromGFF(
            "shu",
            get_sample_data("mbf_genomics/test_with_name.gff3"),
            get_genome_chr_length(genome=get_genome()),
            filter_function=lambda entry: entry["source"]
            == b"Regions_of_sig_enrichment",
            comment_char="#",
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 4
        assert (a.df["chr"] == ["2", "3", "3", "3"]).all()
        assert (a.df["start"] == [662_743, 53968, 58681, 68187]).all()
        assert (
            a.df["score"]
            == [
                b"-0.583140831299777",
                b"0.41667200444801",
                b"1.20483507668127",
                b"0.742024408736128",
            ]
        ).all()
        assert (a.df["name"] == ["reg5", "reg1", "reg2", "reg3"]).all()

    def test_gff_below_zero(self):
        genome = get_genome()
        b = regions.GenomicRegions_FromGFF(
            "sha",
            get_sample_data("mbf_genomics/test_below_zero.gff3"),
            get_genome_chr_length(genome=genome),
            filter_function=lambda entry: entry["source"]
            == b"Regions_of_sig_enrichment",
            comment_char="#",
            fix_negative_coordinates=True,
        )

        if ppg.inside_ppg():
            a = regions.GenomicRegions_FromGFF(
                "shu",
                get_sample_data("mbf_genomics/test_below_zero.gff3"),
                get_genome_chr_length(genome=genome),
                filter_function=lambda entry: entry["source"]
                == b"Regions_of_sig_enrichment",
                comment_char="#",
            )
            force_load(a.load())
            force_load(b.load())
            jobA = a.load()
            b.load()
            with pytest.raises(ppg.RuntimeError):
                ppg.run_pipegraph()
                assert isinstance(jobA.exception, "ValueError")
                assert len(b.df)
        else:
            with pytest.raises(ValueError):
                a = regions.GenomicRegions_FromGFF(
                    "shu",
                    get_sample_data("mbf_genomics/test_below_zero.gff3"),
                    get_genome_chr_length(genome=genome),
                    filter_function=lambda entry: entry["source"]
                    == b"Regions_of_sig_enrichment",
                    comment_char="#",
                )
            b.load()
        assert (b.df["start"] >= 0).all()

    def test_gff_chromosome_mangler(self):
        a = regions.GenomicRegions_FromGFF(
            "shu",
            get_sample_data("mbf_genomics/test.gff3"),
            get_genome_chr_length(genome=get_genome()),
            filter_function=lambda entry: entry["source"]
            == b"Regions_of_sig_enrichment",
            comment_char="#",
            chromosome_mangler=lambda chr: str(int(chr) + 1),
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 4
        assert (a.df["chr"] == ["3", "4", "4", "4"]).all()
        assert (a.df["start"] == [662_743, 53968, 58681, 68187]).all()
        assert (
            a.df["score"]
            == [
                b"-0.583140831299777",
                b"0.41667200444801",
                b"1.20483507668127",
                b"0.742024408736128",
            ]
        ).all()

    def test_bed(self):
        a = regions.GenomicRegions_FromBed(
            "shu",
            get_sample_data("mbf_genomics/test.bed"),
            get_genome_chr_length(genome=get_genome()),
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 6
        assert (a.df["chr"] == ["2", "2", "2", "3", "3", "3"]).all()
        assert (
            a.df["start"] == [356_591, 662_743, 1_842_875, 53968, 58681, 68187]
        ).all()
        should = [
            -0.610_793_694_004_74,
            -0.583_140_831_299_777,
            numpy.nan,
            0.416_672_004_448_01,
            1.204_835_076_681_27,
            0.742_024_408_736_128,
        ]
        assert (a.df["name"] == ["four", "five", "six", "one", "two", "three"]).all()
        for ii, val in enumerate(a.df["score"]):
            if numpy.isnan(should[ii]):
                assert numpy.isnan(val)
            else:
                assert round(abs(val - should[ii]), 7) == 0

    def test_bigbed(self):
        a = regions.GenomicRegions_FromBigBed(
            "shu",
            get_sample_data("mbf_genomics/test.bb"),
            get_genome_chr_length(genome=get_genome()),
        )
        force_load(a.load())
        run_pipegraph()
        assert (
            a.df["start"][:5].values == [199_288, 263_869, 264_202, 264_596, 268_492]
        ).all()

    def test_bed_without_score(self):
        a = regions.GenomicRegions_FromBed(
            "shu",
            get_sample_data("mbf_genomics/test_noscore.bed"),
            get_genome_chr_length(genome=get_genome()),
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 6
        assert (a.df["chr"] == ["2", "2", "2", "3", "3", "3"]).all()
        assert (
            a.df["start"] == [356_591, 662_743, 1_842_875, 53968, 58681, 68187]
        ).all()
        assert not "score" in a.df.columns
        assert (a.df["name"] == ["four", "five", "six", "one", "two", "three"]).all()

    def test_bed_constant_name(self):
        a = regions.GenomicRegions_FromBed(
            "shu",
            get_sample_data("mbf_genomics/test_constant_name.bed"),
            get_genome_chr_length(genome=get_genome()),
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 6
        assert (a.df["chr"] == ["2", "2", "2", "3", "3", "3"]).all()
        assert (
            a.df["start"] == [356_591, 662_743, 1_842_875, 53968, 58681, 68187]
        ).all()
        assert not "name" in a.df.columns

    def test_empty_bed(self):
        with RaisesDirectOrInsidePipegraph(ValueError):
            a = regions.GenomicRegions_FromBed(
                "shu",
                get_sample_data("mbf_genomics/test_empty.bed"),
                get_genome_chr_length(genome=get_genome()),
            )
            force_load(a.load())

    def test_wig(self):
        a = regions.GenomicRegions_FromWig(
            "shu",
            get_sample_data("mbf_genomics/test.wig"),
            get_genome_chr_length(genome=get_genome()),
            enlarge_5prime=2,
            enlarge_3prime=1,
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 5
        assert (a.df["chr"] == ["3", "3", "3", "3", "3"]).all()
        assert (a.df["start"] + 2 == [2, 41, 81, 120, 158]).all()
        assert (a.df["stop"] - 1 == [37, 76, 116, 155, 193]).all()
        assert (
            numpy.abs(
                a.df["score"]
                - [
                    0.150_787_654_060_043,
                    0.158_118_185_847_438,
                    0.165_947_064_796_052,
                    0.175_360_778_350_916,
                    0.191_375_614_608_046,
                ]
            )
            <= 0.0001
        ).any()

    def test_partec(self):
        a = regions.GenomicRegions_FromPartec(
            "shu",
            get_sample_data("mbf_genomics/test_partec.txt"),
            get_genome_chr_length(genome=get_genome()),
        )
        force_load(a.load())
        run_pipegraph()
        assert len(a.df) == 5
        assert (a.df["chr"] == ["1", "2", "2", "4", "5"]).all()
        assert (
            a.df["start"]
            == [119_019_045, 33_327_521, 216_002_580, 45_933_913, 138_423_628]
        ).all()
        assert (
            a.df["stop"]
            == [119_019_245, 33_327_690, 216_002_792, 45_934_151, 138_423_806]
        ).all()

    def test_binned(self):
        genome = get_genome_chr_length(genome=get_genome())
        with pytest.raises(ValueError):
            regions.GenomicRegions_BinnedGenome(genome, 10000, "1")
        a = regions.GenomicRegions_BinnedGenome(genome, 10000, ["1"])
        b = regions.GenomicRegions_BinnedGenome(genome, 10000, new_name="full_bins")
        force_load(a.load())
        force_load(b.load())
        run_pipegraph()
        assert len(a.df) == 10
        assert len(b.df) == 10 + 20 + 30 + 40 + 50

    def test_windowed(self):
        genome = get_genome_chr_length(genome=get_genome())
        a = regions.GenomicRegions_Windows(genome, "win1", 10000, 0, ["1"])
        b = regions.GenomicRegions_Windows(genome, "win2", 10000, 10000, ["1"])
        c = regions.GenomicRegions_Windows(genome, "win3", 10000, 0)
        force_load(a.load())
        force_load(b.load())
        force_load(c.load())
        run_pipegraph()
        assert len(a.df) == 10
        assert len(b.df) == 5
        assert len(c.df) == 10 + 20 + 30 + 40 + 50


@pytest.mark.usefixtures("no_pipegraph")
class TestOutsideOfPipegraph:
    def test_ignores_second_loading(self):
        genome = get_genome_chr_length(genome=get_genome())
        counter = [0]

        def sample_data():
            counter[0] += 1
            return pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "3", "5"],
                    "start": [10, 100, 1000, 10000, 100_000],
                    "stop": [1010, 110, 1110, 11110, 111_110],
                }
            )

        a = regions.GenomicRegions("shu", sample_data, [], genome, on_overlap="merge")
        assert counter[0] == 1  # load is immediate
        a.load()
        assert counter[0] == 1
        a.load()
        assert counter[0] == 1


class TestTagCountAnnotator:
    def test_simple(self, new_pipegraph):
        from .test_genes_anno_tag_counts import MockBam, MockLane
        from mbf.genomics import genes

        mb = MockBam("1", 10050, mode=0)
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
            chr_lengths={"1": 10050},
        )
        a = regions.GenomicRegions_BinnedGenome(genome, 1000, ["1"])
        mb_lane = MockLane(mb, genome)
        anno = genes.anno_tag_counts.GRUnstrandedRust(mb_lane)
        a.add_annotator(anno)
        force_load(a.annotate)
        new_pipegraph.run()
        assert (a.df[anno.columns[0]] == [1000] * 10 + [50]).all()

    def test_both_strands(self, new_pipegraph):
        from .test_genes_anno_tag_counts import MockBam, MockLane
        from mbf.genomics import genes

        mb = MockBam("1", 10050, mode=1)
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
            chr_lengths={"1": 10050},
        )
        a = regions.GenomicRegions_BinnedGenome(genome, 1000, ["1"])
        mb_lane = MockLane(mb, genome)
        anno = genes.anno_tag_counts.GRUnstrandedRust(mb_lane)
        a.add_annotator(anno)
        force_load(a.annotate)
        new_pipegraph.run()
        assert (a.df[anno.columns[0]] == [2000] * 10 + [50 * 2]).all()

    def test_both_strands_stranded(self, new_pipegraph):
        from .test_genes_anno_tag_counts import MockBam, MockLane
        from mbf.genomics import genes

        mb = MockBam("1", 10050, mode=1)
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
            chr_lengths={"1": 10050},
        )
        a = regions.GenomicRegions_BinnedGenome(genome, 1000, ["1"])
        b = regions.GenomicRegions(
            "B", lambda: a.df.assign(strand=1), [a.load()], genome
        )
        mb_lane = MockLane(mb, genome)
        anno = genes.anno_tag_counts.GRStrandedRust(mb_lane)
        b.add_annotator(anno)
        force_load(b.annotate)
        new_pipegraph.run()
        assert (b.df[anno.columns[0]] == [1000] * 10 + [50 * 1]).all()
