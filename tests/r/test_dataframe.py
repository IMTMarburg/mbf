import pandas as pd
from pandas.testing import assert_frame_equal
from mbf.r import convert_dataframe_to_r, convert_dataframe_from_r
import rpy2.robjects as ro
import pytest


class TestDataFrameRoundTripping:
    def test_simple_from(self):
        rdf = ro.r("""data.frame(a=c(0.1, 0.2, 0.3))""")
        pdf = convert_dataframe_from_r(rdf)
        should = pd.DataFrame({"a": [0.1, 0.2, 0.3], "idx": ["1", "2", "3"]}).set_index(
            "idx"
        )
        should.index.name = None
        assert_frame_equal(pdf, should)

    def test_simple(self):
        df = pd.DataFrame({"a": [0.1, 0.2, 0.3]})
        rdf = convert_dataframe_to_r(df)
        pdf = convert_dataframe_from_r(rdf)
        df.index = [str(x) for x in df.index]
        assert_frame_equal(df, pdf)

    def test_index_gets_dropped(self):
        df = pd.DataFrame({"a": [0.1, 0.2, 0.3], "idx": ["a", "b", "c"]}).set_index(
            "idx"
        )
        rdf = convert_dataframe_to_r(df)
        pdf = convert_dataframe_from_r(rdf)
        df.index.name = None
        assert_frame_equal(df, pdf)

    def test_categorical(self):
        df = pd.DataFrame({"a": pd.Categorical(["a", "b", "c", "a"])})
        rdf = convert_dataframe_to_r(df)
        print(rdf)
        pdf = convert_dataframe_from_r(rdf)
        print(pdf)
        df.index = [str(x) for x in df.index]
        assert_frame_equal(df, pdf)

    def test_categorical_ordered(self):
        df = pd.DataFrame({"a": pd.Categorical(["a", "b", "c", "a"], ordered=True)})
        rdf = convert_dataframe_to_r(df)
        pdf = convert_dataframe_from_r(rdf)
        df.index = [str(x) for x in df.index]
        assert_frame_equal(df, pdf)

    def test_raise_on_funny(self):
        df = pd.DataFrame({"aaaa": [(1, 2), (3, 4)]})
        with pytest.raises(ValueError) as e:
            convert_dataframe_to_r(df)
        assert "aaaa" in str(e.value)
