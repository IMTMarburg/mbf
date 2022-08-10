import pytest
from pathlib import Path
import pypipegraph as ppg
from mbf.qualitycontrol import (
    register_qc,
    prune_qc,
    get_qc_jobs,
    QCCollectingJob,
    disable_qc,
    qc_disabled,
)

from mbf.qualitycontrol.testing import assert_image_equal


class TestRegistration:
    def test_registration_and_pruning(self, new_pipegraph):
        with pytest.raises(TypeError):
            register_qc("shu")
        jobA = ppg.FileGeneratingJob("a", lambda: Path("a").write_text("hello"))
        register_qc(jobA)
        print(list(get_qc_jobs()))
        assert jobA in list(get_qc_jobs())
        assert not jobA._pruned
        jobc = register_qc(
            ppg.FileGeneratingJob("c", lambda: Path("b").write_text("hello"))
        )

        def check_prune(job):
            return job.job_id.lower()[-1] == "c"

        prune_qc(check_prune)
        assert jobc in list(get_qc_jobs())
        assert not jobc._pruned
        jobB = register_qc(
            ppg.FileGeneratingJob("b", lambda: Path("b").write_text("hello"))
        )
        assert jobB in list(get_qc_jobs())
        assert jobB._pruned
        jobC = register_qc(
            ppg.FileGeneratingJob("C", lambda: Path("b").write_text("hello"))
        )
        assert not jobC._pruned
        assert len(list(get_qc_jobs())) == 4
        prune_qc()
        assert jobA._pruned
        assert jobB._pruned
        assert jobc._pruned
        assert jobC._pruned
        for j in get_qc_jobs():
            assert j._pruned

    def test_pruning_plotjob(self, new_pipegraph):
        jobA = register_qc(ppg.PlotJob("c.png", lambda: None, lambda: None))
        assert not jobA._pruned
        prune_qc()
        assert jobA._pruned
        assert jobA.cache_job._pruned
        assert jobA.table_job._pruned

    def test_collecting_qc_job(self, new_pipegraph):
        def output(output_filename, elements):
            Path(output_filename).write_text("\n".join((e.name for e in elements)))

        class Collectible:
            def __init__(self, name):
                self.name = name

        job = QCCollectingJob("output", output)
        job.add(Collectible("hello"))
        job2 = QCCollectingJob("output", output)
        assert job2 is job
        job.add(Collectible("world"))
        ppg.run_pipegraph()
        assert Path("output").read_text() == "hello\nworld"

    def test_qc_pruning_is_bound_to_ppg(self, new_pipegraph):
        assert not qc_disabled()
        disable_qc()
        assert qc_disabled()
        new_pipegraph.new_pipegraph()
        assert not qc_disabled()

    def test_no_pipegraph_means_disabled(self, no_pipegraph):
        assert ppg.util.global_pipegraph is None
        assert qc_disabled()

    def test_both_fixtures(self, both_ppg_and_no_ppg_no_qc):
        if not ppg.inside_ppg():
            assert qc_disabled()
        else:
            assert qc_disabled()

    def test_new_pipegraph_no_qc(self, new_pipegraph_no_qc):
        assert qc_disabled()


class TestAssertImagesEqualInside:
    def test_assert_images_equal_inside_class(self):
        assert_image_equal(
            Path(__file__).parent
            / "base_images"
            / "test_qc"
            / "_"
            / "test_assert_images_equal.png"
        )
        with pytest.raises(ValueError):
            assert_image_equal(
                Path(__file__).parent
                / "base_images"
                / "test_qc"
                / "_"
                / "test_assert_images_equal.png",
                "_b",
            )


def test_assert_images_equal():
    assert_image_equal(
        Path(__file__).parent
        / "base_images"
        / "test_qc"
        / "_"
        / "test_assert_images_equal.png"
    )
    with pytest.raises(ValueError) as e:  # here the baseline image does not exist
        assert_image_equal(
            Path(__file__).parent
            / "base_images"
            / "test_qc"
            / "_"
            / "test_assert_images_equal.png",
            "_b",
        )
    # should_path overwrites suffix
    assert_image_equal(
        Path(__file__).parent
        / "base_images"
        / "test_qc"
        / "_"
        / "test_assert_images_equal.png",
        "_b",
        should_path=Path(__file__).parent
        / "base_images"
        / "test_qc"
        / "_"
        / "test_assert_images_equal.png",
    )

    assert "Base_line image not found" in str(e.value)
    with pytest.raises(ValueError) as e:  # here it is different
        assert_image_equal(
            Path(__file__).parent
            / "base_images"
            / "test_qc"
            / "_"
            / "test_assert_images_equal.png",
            suffix="_c",
        )
    assert "Image files did not match" in str(e.value)

    with pytest.raises(IOError) as e:
        assert_image_equal("does not exist")
    assert "not created" in str(e.value)

    # with pytest.raises(ValueError) as e: #here it is different
    with pytest.raises(ValueError) as e:  # here it is different
        assert_image_equal(
            Path(__file__).parent
            / "base_images"
            / "test_qc"
            / "_"
            / "test_assert_images_equal.png",
            suffix="_d",
        )
    assert "do not match expected size" in str(e.value)
