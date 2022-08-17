import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
import pypipegraph as ppg
from mbf.genomics.annotator import Constant, Annotator
from mbf.genomics.util import (
    read_pandas,
    freeze,
    parse_a_or_c_to_column,
    parse_a_or_c_to_anno,
    parse_a_or_c_to_plot_name,
    find_annos_from_column,
)


def test_read_pandas_csv_in_xls(new_pipegraph):
    df = pd.DataFrame({"x": ["a", "b", "c"], "y": [1, 2.5, 3]})
    df.to_excel("shu.xlsx", index=False)
    assert_frame_equal(df, read_pandas("shu.xlsx"))
    df.to_csv("shu.xls", sep="\t", index=False)
    assert_frame_equal(df, read_pandas("shu.xls"))
    df.to_csv("shu.tsv", sep="\t", index=False)
    assert_frame_equal(df, read_pandas("shu.tsv"))
    df.to_csv("shu.csv", index=False)
    assert_frame_equal(df, read_pandas("shu.csv"))
    df.to_csv("shu.something", index=False)
    with pytest.raises(ValueError):
        read_pandas("shu.something")


def test_freeze():
    a = {"a": [1, 2, 3], "b": {"c": set([2, 3, 5])}}
    with pytest.raises(TypeError):
        hash(a)
    assert hash(freeze(a))
    if not hasattr(ppg, "is_ppg2"):
        assert freeze(a) == freeze(freeze(a))

        class Nohash:
            def __hash__(self):
                return NotImplemented

        with pytest.raises(TypeError):
            freeze(Nohash())


class PolyConstant(Annotator):
    def __init__(self, column_names, values, plot_name=None):
        self.columns = column_names
        self.value = values
        if plot_name is not None:
            self.plot_name = plot_name

    def calc(self, df):
        return pd.DataFrame(
            {k: v for (k, v) in zip(self.columns, self.value)}, index=df.index
        )


class TestAnnotatorParsing:
    def test_to_column(self):
        assert parse_a_or_c_to_column("hello") == "hello"
        assert parse_a_or_c_to_column(Constant("shu", 5)) == "shu"
        assert parse_a_or_c_to_column(PolyConstant(["shu", "sha"], [5, 10])) == "shu"
        assert (
            parse_a_or_c_to_column((PolyConstant(["shu", "sha"], [5, 10]), 1)) == "sha"
        )
        assert (
            parse_a_or_c_to_column((PolyConstant(["shu", "sha"], [5, 10]), "sha"))
            == "sha"
        )
        assert parse_a_or_c_to_column((None, "shi")) == "shi"
        with pytest.raises(KeyError):
            parse_a_or_c_to_column((PolyConstant(["shu", "sha"], [5, 10]), "shi"))
        with pytest.raises(IndexError):
            parse_a_or_c_to_column((PolyConstant(["shu", "sha"], [5, 10]), 5))

        with pytest.raises(ValueError):
            parse_a_or_c_to_column(5)
        with pytest.raises(ValueError):
            parse_a_or_c_to_column((Constant("shu", 5), "shu", 3))

    def test_to_anno(self):
        assert parse_a_or_c_to_anno("hello") is None
        assert parse_a_or_c_to_anno(Constant("shu", 5)) == Constant("shu", 5)
        assert parse_a_or_c_to_anno(
            PolyConstant(["shu", "sha"], [5, 10])
        ) == PolyConstant(["shu", "sha"], [5, 10])
        assert parse_a_or_c_to_anno(
            (PolyConstant(["shu", "sha"], [5, 10]), 1)
        ) == PolyConstant(["shu", "sha"], [5, 10])
        assert parse_a_or_c_to_anno(
            (PolyConstant(["shu", "sha"], [5, 10]), "sha")
        ) == PolyConstant(["shu", "sha"], [5, 10])
        with pytest.raises(KeyError):
            parse_a_or_c_to_anno((PolyConstant(["shu", "sha"], [5, 10]), "shi"))
        with pytest.raises(IndexError):
            parse_a_or_c_to_anno((PolyConstant(["shu", "sha"], [5, 10]), 5))

        with pytest.raises(ValueError):
            parse_a_or_c_to_anno(5)
        with pytest.raises(ValueError):
            parse_a_or_c_to_anno((Constant("shu", 5), "shu", 3))

    def test_to_plot_name(self):
        assert parse_a_or_c_to_plot_name("hello") == "hello"
        assert parse_a_or_c_to_plot_name(Constant("shu", 5)) == "shu"
        assert parse_a_or_c_to_plot_name(PolyConstant(["shu", "sha"], [5, 10])) == "shu"
        assert (
            parse_a_or_c_to_plot_name((PolyConstant(["shu", "sha"], [5, 10]), 1))
            == "sha"
        )
        assert (
            parse_a_or_c_to_plot_name((PolyConstant(["shu", "sha"], [5, 10]), "sha"))
            == "sha"
        )
        with pytest.raises(KeyError):
            parse_a_or_c_to_plot_name((PolyConstant(["shu", "sha"], [5, 10]), "shi"))
        with pytest.raises(IndexError):
            parse_a_or_c_to_plot_name((PolyConstant(["shu", "sha"], [5, 10]), 5))

        with pytest.raises(ValueError):
            parse_a_or_c_to_plot_name(5)
        with pytest.raises(ValueError):
            parse_a_or_c_to_plot_name((Constant("shu", 5), "shu", 3))

        assert (
            parse_a_or_c_to_plot_name(PolyConstant(["shu", "sha"], [5, 10], "hello"))
            == "hello"
        )
        assert (
            parse_a_or_c_to_plot_name(
                (PolyConstant(["shu", "sha"], [5, 10], "hello"), "sha")
            )
            == "hello"
        )
        assert (
            parse_a_or_c_to_plot_name(
                (PolyConstant(["shu", "sha"], [5, 10], "hello"), 1)
            )
            == "hello"
        )

    def test_find_annos_from_column(self, both_ppg_and_no_ppg_no_qc, clear_annotators):
        a = Constant("shu", 5)
        assert find_annos_from_column("shu") == [a]
        assert find_annos_from_column("shu")[0] is a
        with pytest.raises(KeyError):
            find_annos_from_column("nosuchcolumn")

        b = PolyConstant(["shu"], [10])
        assert find_annos_from_column("shu") == [a, b]

        if ppg.inside_ppg():
            both_ppg_and_no_ppg_no_qc.new_pipegraph()
            with pytest.raises(KeyError):
                find_annos_from_column("shu")
