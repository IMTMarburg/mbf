import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
import pypipegraph as ppg
from pathlib import Path
from mbf_genomics import DelayedDataFrame
from mbf_genomics.annotator import Annotator
from .shared import force_load_ddf as force_load


def DummyAnnotatable(name):
    return DelayedDataFrame(
        name,
        lambda: pd.DataFrame(
            {
                "a": ["a", "b", "c", "d"],
                "b": [1, 2, 3, 4],
                "c": [200.1, 100.2, 400.3, 300.4],
            }
        ),
    )


class SequenceAnnotator(Annotator):
    columns = ["sequence"]

    def calc(self, df):
        return pd.DataFrame({self.columns[0]: range(0, len(df))})


class SequenceAnnotatorDuo(Annotator):
    columns = ["sequenceDuo", "rev_sequenceDuo"]

    def calc(self, df):
        return pd.DataFrame(
            {self.columns[0]: range(0, len(df)), self.columns[1]: range(len(df), 0, -1)}
        )


class SequenceAnnotatorDuoCollision(Annotator):
    columns = ["shu", "rev_sequenceDuo"]

    def calc(self, df):
        return pd.DataFrame(
            {self.columns[0]: range(0, len(df)), self.columns[1]: range(len(df), 0, -1)}
        )


class FixedAnnotator(Annotator):
    def __init__(self, column_name, values):
        self.columns = [column_name]
        self.values = values

    def deps(self, ddf):
        return ppg.ParameterInvariant(
            ddf.name + "_" + self.columns[0], str(self.values)
        )

    def calc(self, df):
        op = open("dummy.txt", "ab")
        op.write(b"A")
        op.close()
        return pd.DataFrame({self.columns[0]: self.values[: len(df)]})


class FixedAnnotator2(Annotator):  # used for conflict of annotator class tests
    def __init__(self, column_name, values):
        self.columns = [column_name]
        self.values = values

    def deps(self, ddf):
        return ppg.ParameterInvariant(
            ddf.name + "_" + self.column_name, str(self.values)
        )

    def annotate(self, annotat):
        op = open("dummy.txt", "ab")
        op.write(b"A")
        op.close()
        return pd.DataFrame({self.columns[0]: self.values[: len(annotat)]})


class BrokenAnnoDoesntCallAnnotatorInit(Annotator):
    columns = ["shu"]

    def calc(self, df):
        return pd.DataFrame({self.column_name: range(0, len(df))})


class FakeAnnotator(object):

    columns = ["shu"]

    def calc(self, df):
        return pd.DataFrame({self.columns[0]: range(0, len(df))})


@pytest.mark.usefixtures("new_pipegraph")
class Test_FromOldGenomics:
    def test_add_annotator_takes_only_annotators(self):
        a = DummyAnnotatable("A")
        with pytest.raises(TypeError):
            a += 123

    def test_non_anno_raises(self):
        a = DummyAnnotatable("A")
        with pytest.raises(TypeError):
            a += FakeAnnotator()

    def test_one_column_annotator(self, new_pipegraph):
        # new_pipegraph.new(log_level=6)
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        a.add_annotator(anno)
        force_load(a)
        ppg.run_pipegraph()
        assert (a.df["sequence"] == [0, 1, 2, 3]).all()

    def test_two_column_annotator(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotatorDuo()
        a.add_annotator(anno)
        force_load(a)
        ppg.run_pipegraph()
        assert (a.df["sequenceDuo"] == [0, 1, 2, 3]).all()
        assert (a.df["rev_sequenceDuo"] == [4, 3, 2, 1]).all()

    def test_two_differenct_annotators_with_identical_column_names_raise_on_adding(
        self,
    ):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotatorDuo()
        a.add_annotator(anno)
        anno2 = SequenceAnnotatorDuoCollision()
        a.add_annotator(anno2)
        force_load(a)
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()

    def test_annotator_copying_on_filter(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        force_load(even)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["b"] == [2, 4]).all()
        assert (even.df["sequence"] == [1, 3]).all()

    def test_annotator_copying_on_filter_two_deep(self, new_pipegraph, job_trace_log):
        if hasattr(ppg, "is_ppg2"):
            new_pipegraph.new_pipegraph(log_level=6)

            a = DummyAnnotatable("A")
            anno = SequenceAnnotator()
            even = a.filter("event", lambda df: df["b"] % 2 == 0)
            force_load(even)
            second = even.filter("event2", lambda df: df["b"] == 4)
            a.add_annotator(anno)
            force_load(second)
            # second.write('shu.tsv')

            ppg.run_pipegraph()
            assert (second.df["b"] == [4]).all()
            assert (second.df["sequence"] == [3]).all()

    def test_annotator_copying_on_filter_with_anno(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        even = a.filter("event", lambda df: df["sequence"] % 2 == 0, annotators=[anno])
        force_load(even)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["b"] == [1, 3]).all()
        assert (even.df["sequence"] == [0, 2]).all()

    def test_no_anno_data_copying_if_no_annotate_dependency(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)

        def write():
            op = open("dummy.txt", "wb")
            op.write(b"SHU")
            op.close()

        ppg.FileGeneratingJob("dummy.txt", write).depends_on(even.load())
        ppg.run_pipegraph()
        assert (even.df["b"] == [2, 4]).all()
        assert "sequence" not in even.df.columns

    def test_anno_data_copying_if_add_annotator_dependency(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)

        def wf():
            op = open("dummy.txt", "wb")
            op.write(b"SHU")
            op.close()

        fg = ppg.FileGeneratingJob("dummy.txt", wf)
        even.add_annotator(anno)
        fg.depends_on(even.add_annotator(anno))
        ppg.run_pipegraph()
        assert (even.df["b"] == [2, 4]).all()
        assert (even.df["sequence"] == [1, 3]).all()

    def test_annotator_copying_on_sort_and_top(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        a.add_annotator(anno)
        even = a.filter(
            "event", lambda df: df.sort_values("b", ascending=False)[:2].index
        )
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["b"] == [4, 3]).all()
        assert (even.df["sequence"] == [3, 2]).all()

    def test_annotator_just_added_to_child(self):
        a = DummyAnnotatable("A")
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        anno = SequenceAnnotator()
        even.add_annotator(anno)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["b"] == [2, 4]).all()
        # after all, we add it anew.
        assert (even.df["sequence"] == [0, 1]).all()
        assert "sequence" not in a.df.columns

    def test_annotator_first_added_to_parent_then_to_child(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        even.add_annotator(anno)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["b"] == [2, 4]).all()
        assert (even.df["sequence"] == [1, 3]).all()
        assert (a.df["sequence"] == [0, 1, 2, 3]).all()

    def test_annotator_first_added_to_parent_then_to_second_child(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0).filter(
            "shu", lambda df: df["b"] == 2
        )
        even.add_annotator(anno)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["b"] == [2]).all()
        assert (even.df["sequence"] == [1]).all()
        assert (a.df["sequence"] == [0, 1, 2, 3]).all()

    def test_annotator_first_added_to_child_then_to_parent(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        even.add_annotator(anno)
        force_load(even)
        a.add_annotator(anno)
        force_load(a)
        ppg.run_pipegraph()
        assert "sequence" in a.df.columns
        assert "sequence" in even.df.columns

    def test_annotator_added_after_filtering(self):
        a = DummyAnnotatable("A")
        anno = SequenceAnnotator()
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        a.add_annotator(anno)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["b"] == [2, 4]).all()
        assert (even.df["sequence"] == [1, 3]).all()
        assert (a.df["sequence"] == [0, 1, 2, 3]).all()

    def test_non_hashable_init__args(self):
        a = FixedAnnotator("shu", ["h", "i", "j", "k"])
        b = FixedAnnotator("shu", ["h", "i", "j", "k"])
        assert a is b

    def test_annotator_copying_parent_changed(self, new_pipegraph):
        # first run
        a = DummyAnnotatable("A")
        anno = FixedAnnotator("shu", ("h", "i", "j", "k"))
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["shu"] == ["i", "k"]).all()

        assert Path("dummy.txt").read_text() == "A"  # so it ran once...

        new_pipegraph.new_pipegraph()
        a = DummyAnnotatable("A")
        anno = FixedAnnotator("shu", ("h", "i", "j", "k"))
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["shu"] == ["i", "k"]).all()
        assert Path("dummy.txt").read_text() == "A"  # so it was not rerun

        new_pipegraph.new_pipegraph()
        a = DummyAnnotatable("A")
        anno = FixedAnnotator("shu", ("h", "i", "j", "z"))
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        force_load(even)
        ppg.run_pipegraph()
        assert (even.df["shu"] == ["i", "z"]).all()
        assert Path("dummy.txt").read_text() == "AA"  # so it was rerun

    def test_filter_annotator_copy_nested(self):
        # first run
        a = DummyAnnotatable("A")
        a.write()
        anno = FixedAnnotator("shu", ("h", "i", "j", "k"))
        anno2 = FixedAnnotator("shaw", ("a1", "b2", "c3", "d4"))
        a.add_annotator(anno)
        first = a.filter("first", lambda df: (df["a"] == "b") | (df["a"] == "d"))
        second = first.filter("second", lambda df: ([True, True]))
        third = second.filter("third", lambda df: (df["shu"] == "i"), annotators=[anno])
        fourth = first.filter("fourth", lambda df: ([False, True]))
        second.write()
        fn_4 = fourth.write()[1]
        a.add_annotator(anno2)
        fourth.add_annotator(anno2)
        force_load(first)
        force_load(second)
        force_load(third)
        force_load(fourth)
        ppg.run_pipegraph()
        assert (first.df["shu"] == ["i", "k"]).all()
        assert (first.df["parent_row"] == [1, 3]).all()
        assert (first.df["shaw"] == ["b2", "d4"]).all()
        assert (second.df["shu"] == ["i", "k"]).all()
        assert (second.df["parent_row"] == [1, 3]).all()
        assert (second.df["shaw"] == ["b2", "d4"]).all()
        assert (third.df["shu"] == ["i"]).all()
        assert (third.df["shaw"] == ["b2"]).all()
        assert (third.df["parent_row"] == [1]).all()
        assert (fourth.df["shu"] == ["k"]).all()
        assert (fourth.df["parent_row"] == [3]).all()
        assert (fourth.df["shaw"] == ["d4"]).all()
        df = pd.read_csv(fn_4, sep="\t")
        assert (df["shaw"] == ["d4"]).all()
        assert_frame_equal(
            df, fourth.df.reset_index(drop=True), check_exact=False, rtol=1e-3
        )

    def test_changing_anno_that_filtering_doesnt_care_about_does_not_retrigger_child_rebuild(
        self, new_pipegraph
    ):
        def count():
            op = open("dummyZZ.txt", "ab")
            op.write(b"A")
            op.close()

        fg = ppg.FileGeneratingJob("dummyZZ.txt", count)
        a = DummyAnnotatable("A")
        anno = FixedAnnotator("shu", ("h", "i", "j", "k"))
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        fg.depends_on(even.load())
        ppg.run_pipegraph()
        Path("dummyZZ.txt").read_text() == "A"  # so it ran once...

        new_pipegraph.new_pipegraph()
        fg = ppg.FileGeneratingJob("dummyZZ.txt", count)
        a = DummyAnnotatable("A")
        anno = FixedAnnotator("shu", ("h", "i", "j", "z"))
        a.add_annotator(anno)
        even = a.filter("event", lambda df: df["b"] % 2 == 0)
        fg.depends_on(even.load())
        ppg.run_pipegraph()
        Path("dummyZZ.txt").read_text() == "A"  # so it was not rerun!
        pass

    def test_same_annotor_call_returns_same_object(self):
        anno = FixedAnnotator("shu", ("h", "i", "j", "k"))
        anno2 = FixedAnnotator("shu", ("h", "i", "j", "k"))
        assert anno is anno2

    def test_new_pipeline_invalidates_annotor_cache(self, new_pipegraph):
        anno = FixedAnnotator("shu", ("h", "i", "j", "k"))
        new_pipegraph.new_pipegraph()
        anno2 = FixedAnnotator("shu", ("h", "i", "j", "k"))
        assert anno is not anno2

    def test_raises_on_same_column_name_differing_parameters(self):
        a = DummyAnnotatable("A")
        a += FixedAnnotator("shu", ("h", "i", "j", "k"))
        with pytest.raises(ValueError):
            a += FixedAnnotator("shu", ("h", "i", "j", "h"))

    def test_raises_on_same_column_name_different_annotators(self):
        a = DummyAnnotatable("A")
        a += FixedAnnotator("shu", ("h", "i", "j", "k"))
        with pytest.raises(ValueError):
            a += FixedAnnotator2("shu", ("h", "i", "j", "k"))

    def test_write(self):
        a = DummyAnnotatable("A")
        anno = FixedAnnotator("shu", ("h", "i", "j", "z"))
        a.add_annotator(anno)
        a.write(Path("shu.xls").absolute())
        ppg.run_pipegraph()
        df = pd.read_excel("shu.xls")
        assert_frame_equal(df, a.df, check_dtype=False)  # check_less_precise=2,


@pytest.mark.usefixtures("new_pipegraph")
class TestDynamicAnnotators:
    def test_basic(self):
        class DA(Annotator):
            @property
            def columns(self):
                return ["DA1-A"]

            def deps(self, annotatable):
                return ppg.ParameterInvariant(self.columns[0], "hello")

            def calc(self, df):
                ll = len(df)
                return pd.DataFrame({"DA1-A": [0] * ll})

        a = DummyAnnotatable("A")
        anno = DA()
        a.add_annotator(anno)
        force_load(a)
        ppg.run_pipegraph()
        assert "DA1-A" in a.df.columns
        assert (a.df["DA1-A"] == 0).all()

    def test_multiple_columns(self):
        class DA(Annotator):
            @property
            def columns(self):
                return ["DA2-A", "DA2-B"]

            def deps(self, annotatable):
                return ppg.ParameterInvariant(self.columns[0], "hello")

            def calc(self, df):
                ll = len(df)
                return pd.DataFrame({"DA2-A": [0] * ll, "DA2-B": [1] * ll})

        a = DummyAnnotatable("A")
        anno = DA()
        a.add_annotator(anno)
        force_load(a)
        ppg.run_pipegraph()
        assert "DA2-A" in a.df.columns
        assert (a.df["DA2-A"] == 0).all()
        assert "DA2-B" in a.df.columns
        assert (a.df["DA2-B"] == 1).all()
        assert "DA2-C" not in a.df.columns

    def test_two_differenct_annotators_with_identical_column_names_raise_on_creation(
        self,
    ):
        a = DummyAnnotatable("A")
        columns_called = [False]

        class DA(Annotator):
            def __init__(self, prefix):
                self.prefix = prefix
                self.cache_name = prefix

            @property
            def columns(self):
                raise ValueError()
                columns_called[0] = True
                return ["%s-A" % self.prefix]

            def calc(self, df):
                ll = len(df)
                return pd.DataFrame({"DA1-A": [0] * ll})

        class DA2(Annotator):
            cache_name = "DA2"

            def __init__(self, prefix):
                self.prefix = prefix

            @property
            def columns(self):
                columns_called[0] = True
                return ["%s-A" % self.prefix]

            def annotate(self, df):
                ll = len(df)
                return pd.DataFrame({"DA1-A": [0] * ll})

        a += DA("DA-1")
        d = DA("DA-2")
        a += d  # still ok.
        a += d  # still ok...a
        assert DA("DA-2") is d
        assert columns_called[0] is False
        # with pytest.raises(ppg.RuntimeError):
        ppg.run_pipegraph()

    def test_returning_non_prefix_columns_raises(self):
        class DA(Annotator):
            def __init__(self):
                self.prefix = "DA2-"

            @property
            def columns(self):
                return ["DA2-A", "DA2-B"]

            def calc(self, df):
                ll = len(df)
                return pd.DataFrame({"DB2-A": [0] * ll, "DB2-B": [1] * ll})

        a = DummyAnnotatable("A")
        anno = DA()
        a.add_annotator(anno)
        force_load(a)
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        exceptions = []
        for job in ppg.util.global_pipegraph.jobs.values():
            if job.failed:
                if isinstance(job.exception, Exception):
                    exceptions.append(job.exception)
        assert len(exceptions) == 1
        assert "Annotator declared" in str(exceptions[0])
        assert "DB2-A" in str(exceptions[0])
        assert "DB2-B" in str(exceptions[0])
        assert "DA2-" in str(exceptions[0])

    def test_inheritance(self):
        self._inheritance_test_count = 0

        class DA(Annotator):
            @property
            def columns(self):
                return ["DA2-A"]

            def deps(self, annotatable):
                return ppg.ParameterInvariant(self.columns[0], "hello")

            def calc(self_inner, df):
                ll = len(df)
                res = pd.DataFrame({"DA2-A": [self._inheritance_test_count] * ll})
                self._inheritance_test_count += 1
                return res

        a = DummyAnnotatable("A")
        anno = DA()
        a.add_annotator(anno)
        b = a.filter("B", lambda df: df["b"] > 2)
        b.annotate()
        force_load(a)
        force_load(b)
        ppg.run_pipegraph()
        assert (a.df["DA2-A"] == 0).all()
        assert (b.df["DA2-A"] == 0).all()
        assert len(b.df) == 2

    def test_adding_to_child_first(self):
        self._inheritance_test_count = 0

        class DA(Annotator):
            @property
            def columns(self):
                return ["DB2-A"]

            def deps(self, annotatable):
                return ppg.ParameterInvariant(self.columns[0], "hello")

            def calc(inner_self, df):
                ll = len(df)
                res = pd.DataFrame({"DB2-A": [self._inheritance_test_count] * ll})
                self._inheritance_test_count += 1
                return res

        a = DummyAnnotatable("A")
        anno = DA()
        b = a.filter("B", lambda df: df["b"] > 2)
        b.add_annotator(anno)
        a.add_annotator(anno)
        force_load(b)
        ppg.run_pipegraph()
        assert (b.df["DB2-A"] == 0).all()

    def test_regular_anno_add_identical_column_name_first(self):
        class DA(Annotator):
            @property
            def columns(self):
                return ["DA2-A"]

            def calc(inner_self, df):
                ll = len(df)
                res = pd.DataFrame({"DA2-A": [0] * ll})
                return res

        class RegularAnnotator(Annotator):
            columns = ["DA2-A"]

            def calc(self, df):
                return pd.DataFrame({self.column_name: range(0, len(df))})

        a = DummyAnnotatable("A")
        anno = DA()
        anno2 = RegularAnnotator()
        a.add_annotator(anno2)
        with pytest.raises(ValueError):
            a.add_annotator(anno)

    def test_regular_anno_add_identical_column_name_second(self):
        class DA(Annotator):
            @property
            def columns(self):
                return ["DA2-A"]

            def calc(inner_self, df):
                ll = len(df)
                res = pd.DataFrame({"DA2-A": [0] * ll})
                return res

        class RegularAnnotator(Annotator):
            columns = ["DA2-A"]

            def calc(self, df):
                return pd.DataFrame({self.column_name: range(0, len(df))})

        a = DummyAnnotatable("A")
        anno = DA()
        anno2 = RegularAnnotator()
        a.add_annotator(anno)
        with pytest.raises(ValueError):  # since they share a cache name
            a.add_annotator(anno2)
