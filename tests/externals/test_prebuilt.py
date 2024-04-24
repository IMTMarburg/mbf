import pytest
import pypipegraph as ppg
from pathlib import Path
from mbf.externals import PrebuildManager
from mbf.externals.util import UpstreamChangedError


def counter(filename):
    """Helper for counting invocations in a side-effect file"""
    try:
        res = int(Path(filename).read_text())
    except:  # noqa: E722
        res = 0
    Path(filename).write_text(str(res + 1))
    return str(res)


class TestPrebuilt:
    def test_simple(self, new_pipegraph):
        new_pipegraph.quiet = False
        Path("prebuilt").mkdir()
        mgr = PrebuildManager("prebuilt", "test_host")
        input_files = [Path("one"), Path("two")]
        input_files[0].write_text("hello")
        input_files[1].write_text("world")
        output_files = [Path("outA")]
        count_file = Path("count")
        count_file.write_text("0")

        def calc(output_path):
            counter("calc")
            t = "\n".join([i.read_text() for i in input_files])
            c = int(count_file.read_text())
            (output_path / output_files[0]).write_text(t + str(c))
            count_file.write_text(str(c + 1))

        jobA = mgr.prebuild("dummy", "0.1", input_files, output_files, calc)
        ppg.FileGeneratingJob(
            "shu",
            lambda: Path("shu").write_text(Path(jobA.find_file("outA")).read_text()),
        ).depends_on(jobA)
        ppg.util.global_pipegraph.run()
        # assert Path("prebuilt/test_host/dummy/0.1/outA").read_text() == "hello\nworld0"
        # assert Path("prebuilt/test_host/dummy/0.1/outA.md5sum").exists()
        assert jobA.find_file("outA").read_text() == "hello\nworld0"
        if not hasattr(ppg, "is_ppg2"):
            assert Path("prebuilt/test_host/dummy/0.1/outA.md5sum").exists()
        assert Path("shu").read_text() == "hello\nworld0"
        assert Path("calc").read_text() == "1"

        # no rerunning.
        new_pipegraph = new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild("dummy", "0.1", input_files, output_files, calc)
        ppg.FileGeneratingJob(
            "shu",
            lambda: Path("shu").write_text(Path(jobA.find_file("outA")).read_text()),
        ).depends_on(jobA)
        ppg.util.global_pipegraph.run()
        assert jobA.find_file("outA").read_text() == "hello\nworld0"
        assert Path("shu").read_text() == "hello\nworld0"
        assert Path("calc").read_text() == "1"

        if not hasattr(ppg, "is_ppg2"):
            # no rerunning, getting from second path...
            new_pipegraph = new_pipegraph.new_pipegraph()
            mgr = PrebuildManager("prebuilt", "test_host2")
            assert not Path("prebuilt/test_host2/dummy/0.1/").exists()
            jobA = mgr.prebuild("dummy", "0.1", input_files, output_files, calc)
            assert not Path("prebuilt/test_host2/dummy/0.1/").exists()
            ppg.FileGeneratingJob(
                "shu2",
                lambda: Path("shu2").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            ppg.util.global_pipegraph.run()
            assert not Path("prebuilt/test_host2/dummy/0.1/").exists()
            assert (
                Path("prebuilt/test_host/dummy/0.1/outA").read_text() == "hello\nworld0"
            )
            assert Path("shu2").read_text() == "hello\nworld0"

        if hasattr(ppg, "is_ppg2"):
            new_pipegraph.new_pipegraph()
            mgr = PrebuildManager("prebuilt", "test_host")
            input_files[1].write_text("world!")
            jobA = mgr.prebuild(
                "dummy", "0.1", input_files, output_files, calc, remove_unused=False
            )
            ppg.FileGeneratingJob(
                "shu3",
                lambda: Path("shu3").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            ppg.util.global_pipegraph.run()
            assert Path("shu3").read_text() == "hello\nworld!1"
            assert Path("calc").read_text() == "2"
            assert (
                len(
                    list(
                        Path(mgr.prebuilt_path / mgr.hostname / "dummy" / "done").glob(
                            "*"
                        )
                    )
                )
                == 2
            )  # no cleanup

            new_pipegraph.new_pipegraph()
            mgr = PrebuildManager("prebuilt", "test_host")
            input_files[1].write_text("world!!")
            jobA = mgr.prebuild(
                "dummy",
                "0.1",
                input_files,
                output_files,
                calc,  # remove_unused=True - default
            )
            ppg.FileGeneratingJob(
                "shu3",
                lambda: Path("shu3").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            ppg.util.global_pipegraph.run()
            assert Path("shu3").read_text() == "hello\nworld!!2"
            assert Path("calc").read_text() == "3"  # not a rerun...
            assert (
                len(
                    list(
                        Path(mgr.prebuilt_path / mgr.hostname / "dummy" / "done").glob(
                            "*"
                        )
                    )
                )
                == 1
            )  # cleanup! but only if the job actually runs

        else:
            # changes to the input files, same machine -> explode.
            new_pipegraph.new_pipegraph()
            mgr = PrebuildManager("prebuilt", "test_host")
            input_files[1].write_text("world!")
            jobA = mgr.prebuild("dummy", "0.1", input_files, output_files, calc)
            ppg.FileGeneratingJob(
                "shu3",
                lambda: Path("shu3").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            with pytest.raises(UpstreamChangedError):
                ppg.util.global_pipegraph.run()
            assert not Path("shu3").exists()
            assert (
                Path("prebuilt/test_host/dummy/0.1/outA").read_text() == "hello\nworld0"
            )

        if not hasattr(ppg, "is_ppg2"):
            # changes to the input files, different machine -> explode.
            new_pipegraph.new_pipegraph()
            mgr = PrebuildManager("prebuilt", "test_host2")
            input_files[1].write_text("world!")
            jobA = mgr.prebuild("dummy", "0.1", input_files, output_files, calc)
            ppg.FileGeneratingJob(
                "shu3",
                lambda: Path("shu3").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            with pytest.raises(UpstreamChangedError):
                ppg.util.global_pipegraph.run()
            assert not Path("shu3").exists()
            assert (
                Path("prebuilt/test_host/dummy/0.1/outA").read_text() == "hello\nworld0"
            )
            with pytest.raises(KeyError):
                jobA.find_file("not_found")

        # but new version is ok...
        new_pipegraph.new_pipegraph()
        mgr = PrebuildManager("prebuilt", "test_host")
        jobA = mgr.prebuild("dummy", "0.2", input_files, output_files, calc)
        ppg.FileGeneratingJob(
            "shu3",
            lambda: Path("shu3").write_text(Path(jobA.find_file("outA")).read_text()),
        ).depends_on(jobA)
        if hasattr(ppg, "is_ppg2"):
            before = list(
                Path(mgr.prebuilt_path / mgr.hostname / "dummy" / "by_input").glob("*")
            )
            ppg.util.global_pipegraph.run()
            assert len(before) == 1
            assert jobA.find_file("outA").read_text() == "hello\nworld!!3"
            assert Path("shu3").read_text() == "hello\nworld!!3"
            after = list(
                Path(mgr.prebuilt_path / mgr.hostname / "dummy" / "by_input").glob("*")
            )
            assert len(after) == 1
            # now we have bumped the version - new done_* dir.
            # but it's a symlink. So a) it must not remove the old one
            # and b) it should be a symlink
            new = [x for x in after if x not in before]
            assert len(new) == 1
            assert new[0].is_symlink()

        else:
            ppg.util.global_pipegraph.run()
            assert jobA.find_file("outA").read_text() == "hello\nworld!1"
            assert Path("shu3").read_text() == "hello\nworld!1"

        if hasattr(ppg, "is_ppg2"):
            pass
        else:
            # request same files with minimum_acceptable_version-> no rebuild...
            new_pipegraph.new_pipegraph()
            jobA = mgr.prebuild(
                "dummy",
                "0.3",
                input_files,
                output_files,
                calc,
                minimum_acceptable_version="0.1",
            )
            ppg.FileGeneratingJob(
                "shu4",
                lambda: Path("shu4").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            ppg.util.global_pipegraph.run()
            assert Path("shu4").read_text() == "hello\nworld!1"

            # but with minimum = version, and only older available -> rebuild
            new_pipegraph.new_pipegraph()
            jobA = mgr.prebuild(
                "dummy",
                "0.3",
                input_files,
                output_files,
                calc,
                minimum_acceptable_version="0.3",
            )
            ppg.FileGeneratingJob(
                "shu5",
                lambda: Path("shu5").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            ppg.util.global_pipegraph.run()
            assert Path("shu5").read_text() == "hello\nworld!2"

        # changing the function leads to an exception (ppg1)
        # or a rebuild on ppg2
        new_pipegraph.new_pipegraph()
        mgr = PrebuildManager("prebuilt", "test_host")

        def calc2(output_path):
            t = "\n".join([i.read_text() for i in input_files])
            c = int(count_file.read_text())
            (output_path / output_files[0]).write_text(t + str(c))
            count_file.write_text(str(c + 1) * 2)

        jobA = mgr.prebuild(
            "dummy",
            "0.3",
            input_files,
            output_files,
            calc2,
            minimum_acceptable_version="0.3",
        )
        ppg.FileGeneratingJob(
            "shu5",
            lambda: Path("shu5").write_text(Path(jobA.find_file("outA")).read_text()),
        ).depends_on(jobA)

        if hasattr(ppg, "is_ppg2"):
            before = list(
                Path(mgr.prebuilt_path / mgr.hostname / "dummy" / "by_input").glob("*")
            )
            ppg.run_pipegraph()
            after = list(
                Path(mgr.prebuilt_path / mgr.hostname / "dummy" / "by_input").glob("*")
            )
            assert len(before) == 1
            assert len(after) == 1
            assert before != after
        else:
            with pytest.raises(UpstreamChangedError):
                ppg.util.global_pipegraph.run()
                raise ValueError(
                    [x.job_id for x in ppg.util.global_pipegraph.jobs.values()]
                )

        # this is also true if it was previously build on another machine
        # only ppg1 does the steal-it-from-another-machine magic
        if hasattr(ppg, "is_ppg2"):
            pass
        else:
            new_pipegraph.new_pipegraph()
            mgr = PrebuildManager("prebuilt", "test_host2")

            def calc2(output_path):
                t = "\n".join([i.read_text() for i in input_files])
                c = int(count_file.read_text())
                (output_path / output_files[0]).write_text(t + str(c))
                count_file.write_text(str(c + 1) * 2)

            jobA = mgr.prebuild(
                "dummy",
                "0.3",
                input_files,
                output_files,
                calc2,
                minimum_acceptable_version="0.3",
            )
            ppg.FileGeneratingJob(
                "shu5",
                lambda: Path("shu5").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            with pytest.raises(UpstreamChangedError):
                ppg.util.global_pipegraph.run()

            # but going back to the original -> ok
            new_pipegraph.new_pipegraph()
            mgr = PrebuildManager("prebuilt", "test_host2")

            def calc3(output_path):
                t = "\n".join([i.read_text() for i in input_files])
                c = int(count_file.read_text())
                (output_path / output_files[0]).write_text(t + str(c))
                count_file.write_text(str(c + 1))

            jobA = mgr.prebuild(
                "dummy",
                "0.3",
                input_files,
                output_files,
                calc,
                minimum_acceptable_version="0.3",
            )
            ppg.FileGeneratingJob(
                "shu5",
                lambda: Path("shu5").write_text(
                    Path(jobA.find_file("outA")).read_text()
                ),
            ).depends_on(jobA)
            ppg.util.global_pipegraph.run()
            assert (
                Path("shu5").read_text() == "hello\nworld!2"
            )  # not rerun, neither the function nor the input files changed

    def test_chained(self, new_pipegraph):
        Path("prebuilt").mkdir()
        mgr = PrebuildManager("prebuilt", "test_host")

        def calc_a(output_path):
            (output_path / "A").write_text("hello")

        jobA = mgr.prebuild("partA", "0.1", [], "A", calc_a)

        def calc_b(output_path):
            (output_path / "B").write_text(jobA.find_file("A").read_text() + " world")

        if hasattr(ppg, "is_ppg2"):
            jobB = mgr.prebuild(
                "partB", "0.1", [], "B", calc_b
            )  # that's not an input file, just add the dependency.
        else:
            jobB = mgr.prebuild("partB", "0.1", [jobA.find_file("A")], "B", calc_b)
        jobB.depends_on(jobA)
        ppg.util.global_pipegraph.run()
        assert jobB.find_file("B").read_text() == "hello world"

    def test_chained2(self, new_pipegraph):
        Path("prebuilt").mkdir()
        count_file = Path("count")
        count_file.write_text("0")
        count_file_b = Path("countb")
        count_file_b.write_text("0")
        mgr = PrebuildManager("prebuilt", "test_host")

        def calc_a(output_path):
            c = int(count_file.read_text())
            (output_path / "A").write_text("hello" + str(c))
            count_file.write_text(str(c + 1))

        jobA = mgr.prebuild("partA", "0.1", [], "A", calc_a)
        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "hello0"
        assert count_file.read_text() == "1"

        new_pipegraph.new_pipegraph()

        def calc_b(output_path):
            c = int(count_file_b.read_text())
            count_file_b.write_text(str(c + 1))
            (output_path / "B").write_text(jobA.find_file("A").read_text() + " world")

        jobA = mgr.prebuild("partA", "0.1", [], "A", calc_a)

        if hasattr(ppg, "is_ppg2"):
            jobB = mgr.prebuild("partB", "0.1", [], "B", calc_b)  # no need,
        else:
            jobB = mgr.prebuild("partB", "0.1", [jobA.find_file("A")], "B", calc_b)
        jobB.depends_on(jobA)
        ppg.util.global_pipegraph.run()
        assert jobB.find_file("B").read_text() == "hello0 world"
        assert count_file.read_text() == "1"
        if hasattr(ppg, "is_ppg2"):
            jobA.depends_on_params("shu")
            assert count_file_b.read_text() == "1"  # runs...
            ppg.run_pipegraph()
            assert count_file.read_text() == "2"  # runs...
            assert count_file_b.read_text() == "2"  # runs...
            assert jobA.find_file("A").read_text() == "hello1"
            assert jobB.find_file("B").read_text() == "hello1 world"

    def test_mixed_state_raises(self, new_pipegraph):
        if not hasattr(ppg, "is_ppg2"):
            Path("prebuilt").mkdir()
            mgr = PrebuildManager("prebuilt", "test_host")

            def calc_a(output_path):
                (output_path / "A").write_text("hello")
                (output_path / "B").write_text("hello")

            jobA = mgr.prebuild("partA", "0.1", [], ["A", "B"], calc_a)
            jobA.find_file("A").write_text("something")
            with pytest.raises(ValueError):
                ppg.util.global_pipegraph.run()

    def test_prebuilt_job_raises_on_non_iterable(self, new_pipegraph):
        from mbf.externals.prebuild import PrebuildJob

        with pytest.raises(TypeError):
            PrebuildJob(5, lambda: 5, "shu")
        with pytest.raises(TypeError):
            PrebuildJob([5], lambda: 5, "shu")
        if hasattr(ppg, "is_ppg2"):
            expected = TypeError
        else:
            expected = ValueError
        with pytest.raises(expected):
            PrebuildJob([Path("shu").absolute()], lambda: 5, "shu")

    def test_minimal_and_maximal_versions(self, new_pipegraph):
        Path("prebuilt").mkdir()
        count_file = Path("count")
        count_file.write_text("0")
        mgr = PrebuildManager("prebuilt", "test_host")

        def calc_04(output_path):
            (output_path / "A").write_text("0.4")
            c = int(count_file.read_text())
            count_file.write_text(str(c + 1))

        def calc_05(output_path):
            (output_path / "A").write_text("0.5")
            c = int(count_file.read_text())
            count_file.write_text(str(c + 1))

        def calc_06(output_path):
            (output_path / "A").write_text("0.6")
            c = int(count_file.read_text())
            count_file.write_text(str(c + 1))

        def calc_07(output_path):
            (output_path / "A").write_text("0.7")
            c = int(count_file.read_text())
            count_file.write_text(str(c + 1))

        jobA = mgr.prebuild("partA", "0.5", [], "A", calc_05)
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA)
        ppg.util.global_pipegraph.quiet = False
        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.5"
        assert Path("checkme").read_text() == "0.5"
        assert count_file.read_text() == "1"

        # no rerun here
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.5",
            [],
            "A",
            calc_05,
            minimum_acceptable_version="0.3",
            maximum_acceptable_version="0.6",
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.5"
        assert Path("checkme").read_text() == "0.5"
        assert count_file.read_text() == "1"

        # no rerun on I want this exact version
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.5",
            [],
            "A",
            calc_05,
            minimum_acceptable_version="0.5",
            maximum_acceptable_version="0.5",
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.5"
        assert Path("checkme").read_text() == "0.5"
        assert count_file.read_text() == "1"

        # but we don't have this one
        new_pipegraph.new_pipegraph()
        ppg.util.global_pipegraph.quiet = False
        jobA = mgr.prebuild(
            "partA",
            "0.6",
            [],
            "A",
            calc_06,
            minimum_acceptable_version="0.6",
            maximum_acceptable_version=None,
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.6"
        assert Path("checkme").read_text() == "0.6"
        assert count_file.read_text() == "2"

        # again no rerun
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.6",
            [],
            "A",
            calc_06,
            minimum_acceptable_version="0.5",
            maximum_acceptable_version=None,
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.6"
        assert Path("checkme").read_text() == "0.6"
        assert jobA.version == "0.6"  # since 0.6 was build
        assert count_file.read_text() == "2"

        # get an older one
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.4",
            [],
            "A",
            calc_04,
            minimum_acceptable_version=None,
            maximum_acceptable_version="0.4",
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.4"
        assert Path("checkme").read_text() == "0.4"
        assert count_file.read_text() == "3"

        # you want 0.4-.. you get' the 0.6, since the build function did not change.
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.7",
            [],
            "A",
            calc_06,  # no change here...
            minimum_acceptable_version="0.4",
            maximum_acceptable_version=None,
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        assert count_file.read_text() == "3"  # no rerun of the build
        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.6"
        assert Path("checkme").read_text() == "0.6"
        if not hasattr(ppg, "is_ppg2"):
            assert jobA.version == "0.6"
            assert count_file.read_text() == "3"  # no rerun of the build
        else:
            assert jobA.version == "0.7"
            # which is guranteed to be the same output as 0.6, otherwise
            # the SharedMultiFileGeneratingJob would have rebuild
            assert (
                count_file.read_text() == "4"
            )  # rebuild because we changed from calc_04 to calc_06

        # you want 0.4-.. but you changed the build func -> 0.7
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.7",
            [],
            "A",
            calc_07,  # this changed
            minimum_acceptable_version="0.4",
            maximum_acceptable_version=None,
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        # and no rebuild
        assert jobA.find_file("A").read_text() == "0.7"
        assert Path("checkme").read_text() == "0.7"
        assert jobA.version == "0.7"
        if hasattr(ppg, "is_ppg2"):
            assert count_file.read_text() == "5"
        else:
            assert count_file.read_text() == "4"

        # you want 0.5 min, with an 05 build func, you get 0.5 - and a rerun
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.5",
            [],
            "A",
            calc_05,
            minimum_acceptable_version="0.5",
            maximum_acceptable_version=None,
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.5"
        assert Path("checkme").read_text() == "0.5"
        assert jobA.version == "0.5"
        if hasattr(ppg, "is_ppg2"):
            assert count_file.read_text() == "6"
        else:
            assert count_file.read_text() == "5"  # was 4?

        # and at last, we want 0.5
        new_pipegraph.new_pipegraph()
        jobA = mgr.prebuild(
            "partA",
            "0.5",
            [],
            "A",
            calc_05,
            minimum_acceptable_version="0.5",
            maximum_acceptable_version="0.5",
        )
        ppg.FileGeneratingJob(
            "checkme",
            lambda: Path("checkme").write_text(jobA.find_file("A").read_text()),
        ).depends_on(jobA).depends_on_params(jobA.version)

        ppg.util.global_pipegraph.run()
        assert jobA.find_file("A").read_text() == "0.5"
        assert Path("checkme").read_text() == "0.5"
        assert jobA.version == "0.5"
        if hasattr(ppg, "is_ppg2"):
            assert count_file.read_text() == "6"  # was 4?
        else:
            assert count_file.read_text() == "5"  # was 4?

    def test_prebuild_job_raises_on_executing_if_dep_on_anything_not_prebuild(
        self, new_pipegraph
    ):
        if not hasattr(ppg, "is_ppg2"):
            Path("prebuilt").mkdir()
            count_file = Path("count")
            count_file.write_text("0")
            mgr = PrebuildManager("prebuilt", "test_host")

            def calc_05(output_path):
                (output_path / "A").write_text("0.5")
                c = int(count_file.read_text())
                count_file.write_text(str(c + 1))

            jobA = mgr.prebuild("partA", "0.5", [], "A", calc_05)
            with pytest.raises(ppg.JobContractError):
                jobA.depends_on(ppg.FunctionInvariant("shu", lambda: 5))

    def test_depends_on_file(self, new_pipegraph):
        from mbf.externals.prebuild import _PrebuildFileInvariantsExploding

        Path("prebuilt").mkdir()
        count_file = Path("count")
        count_file.write_text("0")
        mgr = PrebuildManager("prebuilt", "test_host")

        def calc_05(output_path):
            (output_path / "A").write_text("0.5")
            c = int(count_file.read_text())
            count_file.write_text(str(c + 1))

        jobA = mgr.prebuild("partA", "0.5", [], "A", calc_05)
        jobA.depends_on_file(count_file)
        if hasattr(ppg, "is_ppg2"):
            for p in jobA.upstreams:
                if hasattr(p, "files") and count_file in p.files:
                    break
            else:
                assert False
        else:
            for p in jobA.prerequisites:
                if isinstance(p, _PrebuildFileInvariantsExploding):
                    if count_file in p.filenames:
                        break
            else:
                assert False

    def test_machine_path(self, new_pipegraph):
        pd = Path("prebuild").absolute()
        pd.mkdir(exist_ok=True)
        mgr = PrebuildManager(pd, "test_host")

        def calc_05(output_path):
            Path(output_path / "A").write_text("0.5")

        job = mgr.prebuild("partA", "0.5", [], "A", calc_05)
        ppg.run_pipegraph()
        if hasattr(ppg, "is_ppg2"):
            assert "/test_host/" in str(job["A"])
            assert job["A"].read_text() == "0.5"
        else:
            assert Path("prebuild/test_host/partA/0.5/A").read_text() == "0.5"


class TestPrebuiltOutsideOfPPG:
    def test_prebuild(self, new_pipegraph):
        def calc_05(output_path):
            (output_path / "A").write_text("0.5")

        Path("prebuilt").mkdir()
        Path("count").write_text("0")

        if hasattr(ppg, "is_ppg2"):
            mgr = PrebuildManager("prebuilt", "test_host")
            jobA = mgr.prebuild("partA", "0.5", [], "A", calc_05)
            jobA.depends_on_file(Path("count"))
            jobA.depends_on_func("shu", lambda: None)

            ppg.run_pipegraph()  # make sure it's there.
        else:
            Path("prebuilt/test_host/partA/0.5").mkdir(parents=True, exist_ok=True)
            Path("prebuilt/test_host/partA/0.5/A").write_text("0.5")

        ppg.util.global_pipegraph = None
        mgr = PrebuildManager("prebuilt", "test_host")
        jobA = mgr.prebuild("partA", "0.5", [], "A", calc_05)
        jobA.depends_on_file(Path("count"))
        jobA.depends_on_func("shu", lambda: None)
        assert jobA.find_file("A")
        with pytest.raises(KeyError):
            jobA.find_file("B")
        assert next(iter(jobA)) is jobA

        with pytest.raises(ValueError):
            mgr.prebuild("partB", "0.5", [], "A", calc_05)
