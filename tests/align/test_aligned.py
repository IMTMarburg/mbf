import pytest
import shutil
from pathlib import Path
import pypipegraph as ppg
import mbf.align
import pysam
from mbf.qualitycontrol.testing import assert_image_equal
from mbf.sampledata import get_sample_data, get_sample_path
from mbf.qualitycontrol import prune_qc, get_qc_jobs


def is_ppg2():
    return hasattr(ppg, "is_ppg2")


class DummyGenome:
    name = "Dummy_genome"

    def job_genes(self):
        return None

    def job_transcripts(self):
        return None

    def download_genome(self):
        return []


@pytest.mark.usefixtures("new_pipegraph_no_qc")
class TestAligned:
    def test_from_existing_bam(self):
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
        bam_job = ppg.FileInvariant(bam_path)
        genome = DummyGenome()
        lane = mbf.align.AlignedSample("test_lane", bam_job, genome, False, "AA123")
        assert lane.name == "test_lane"
        assert lane.load()[0] is bam_job
        assert isinstance(lane.load()[1], ppg.FileInvariant)
        assert lane.genome is genome
        assert not lane.is_paired
        assert lane.vid == "AA123"

        with pytest.raises(ValueError):
            mbf.align.AlignedSample("test_lane", bam_job, genome, False, "AA123")
        lane2 = mbf.align.AlignedSample("test_lane2", bam_job, genome, True, "AA123")
        assert lane2.is_paired

        b = lane.get_bam()
        assert isinstance(b, pysam.Samfile)
        b = lane.get_unique_aligned_bam()
        assert isinstance(b, pysam.Samfile)
        assert lane.get_bam_names()[0] == bam_path
        assert lane.get_bam_names()[1] == bam_path + ".bai"

        assert lane.mapped_reads() == 8
        assert lane.unmapped_reads() == 0
        for job in get_qc_jobs():
            assert job._pruned

    def test_lane_invariants_on_non_accepted_value(self):
        genome = DummyGenome()
        with pytest.raises(ValueError):
            mbf.align.AlignedSample("test_lane", 123, genome, False, "AA123")

    def test_lane_raises_on_multifilegeneratingJobWithTwoBams(self):
        mfg = ppg.MultiFileGeneratingJob(["a.bam", "b.bam"], lambda: 5)
        genome = DummyGenome()
        with pytest.raises(ValueError):
            mbf.align.AlignedSample("test_lane", mfg, genome, False, "AA123")

    def test_lane_raises_on_multifilegeneratingJobWithTwoBais(self):
        mfg = ppg.MultiFileGeneratingJob(["a.bam", "a.bam.bai", "b.bam.bai"], lambda: 5)
        genome = DummyGenome()
        with pytest.raises(ValueError):
            mbf.align.AlignedSample("test_lane", mfg, genome, False, "AA123")

    def test_lane_raises_on_multifilegeneratingJobWithNoBAM(self):
        mfg = ppg.MultiFileGeneratingJob(["a.sam"], lambda: 5)
        genome = DummyGenome()
        with pytest.raises(ValueError):
            mbf.align.AlignedSample("test_lane", mfg, genome, False, "AA123")

    def test_lane_invariants_on_string(self):
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
        genome = DummyGenome()
        lane = mbf.align.AlignedSample("test_lane", bam_path, genome, False, "AA123")
        assert isinstance(lane.load()[0], ppg.FileInvariant)

    def test_missing_index_file(self):
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
        no_index = "noindex.bam"
        shutil.copy(bam_path, no_index)
        genome = DummyGenome()
        lane = mbf.align.AlignedSample("test_lane", no_index, genome, False, "AA123")
        assert isinstance(lane.load()[0], ppg.FileInvariant)
        assert isinstance(lane.load()[1], ppg.FileGeneratingJob)
        assert lane.load()[1].job_id != "noindex.bam.bai"
        assert lane.load()[0] in lane.load()[1].prerequisites
        with pytest.raises(FileNotFoundError):
            lane.mapped_reads()
        ppg.run_pipegraph()
        assert lane.mapped_reads() == 8

    def test_creating_index_for_fg_job(self):
        def gen():
            shutil.copy(get_sample_data(Path("mbf_align/ex2.bam")), "sample.bam")

        ppg.util.global_pipegraph.quiet = False

        job = ppg.FileGeneratingJob("sample.bam", gen)
        genome = DummyGenome()
        lane = mbf.align.AlignedSample("test_lane", job, genome, False, "AA123")
        assert isinstance(lane.load()[1], ppg.FileGeneratingJob)
        assert lane.load()[0] in lane.load()[1].prerequisites
        ppg.run_pipegraph()
        assert Path("sample.bam").exists()
        assert Path("sample.bam.bai").exists()

    def test_creating_index_for_mfg(self):
        def gen():
            shutil.copy(get_sample_data(Path("mbf_align/ex2.bam")), "sample.bam")

        ppg.util.global_pipegraph.quiet = False

        job = ppg.MultiFileGeneratingJob(["sample.bam"], gen)
        genome = DummyGenome()
        lane = mbf.align.AlignedSample("test_lane", job, genome, False, "AA123")
        assert isinstance(lane.load()[1], ppg.FileGeneratingJob)
        assert lane.load()[0] in lane.load()[1].prerequisites
        ppg.run_pipegraph()
        assert Path("sample.bam").exists()
        assert Path("sample.bam.bai").exists()

    def test_subtraction_by_read(self):
        from mbf.sampledata import get_human_22_fake_genome

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
            "AA124",
        )  # index creation is automatic
        lane3 = mbf.align.AlignedSample(
            "test_lane3",
            get_sample_data(Path("mbf_align/chipseq_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        lane3_subset = mbf.align.AlignedSample(
            "test_lane3_subset",
            get_sample_data(Path("mbf_align/chipseq_chr22_subset.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic

        lane_empty = lane.post_process(
            mbf.align.post_process.SubtractOtherLane(lane2), new_name="empty"
        )
        lane_full = lane.post_process(
            mbf.align.post_process.SubtractOtherLane(lane3), new_name="full"
        )
        assert lane_empty.result_dir != lane_full.result_dir
        lane_some = lane3.post_process(
            mbf.align.post_process.SubtractOtherLane(lane3_subset),
            result_dir="results/aligned/shu",
        )
        qc_jobs = [lane_some.post_processor_qc_jobs, lane_full.post_processor_qc_jobs]
        prune_qc(lambda job: job in qc_jobs)
        ppg.run_pipegraph()
        assert Path(lane_empty.get_bam_names()[1]).exists()
        assert Path(lane_full.get_bam_names()[1]).exists()
        assert lane_empty.mapped_reads() == 0
        assert lane_full.mapped_reads() == lane.mapped_reads()
        assert lane.mapped_reads() != 0
        assert (
            lane_some.mapped_reads()
            == lane3.mapped_reads() - lane3_subset.mapped_reads()
        )
        assert lane3_subset.mapped_reads()  # make sure there was something to subtract
        assert "shu" in lane_some.get_bam_names()[0]
        assert_image_equal(qc_jobs[0].filenames[0], "_result_dir")
        assert_image_equal(qc_jobs[0].filenames[0])

    def test_to_fastq(self):
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
        bam_job = ppg.FileInvariant(bam_path)
        genome = DummyGenome()
        lane = mbf.align.AlignedSample("test_lane", bam_job, genome, False, "AA123")
        fastq_path = "out.fastq"
        lane.to_fastq(fastq_path)
        ppg.run_pipegraph()
        assert Path(fastq_path).exists()
        assert (
            Path(fastq_path).read_text()
            == """@read_28833_29006_6945
AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG
+
<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<
@read_28701_28881_323b
TGCAAGGCCGCATCGGCCAAGGCCAAGATATAGGT
+
<<<<7<<<<<<<<<<<<;6<<<:;7<<<<;<<<<<
@read_28701_28881_323c
TGCAAGGCCGCATCGGCCAAGGCCAAGATATAGGT
+
<<<<7<<<<<<<<<<<<;6<<<:;7<<<<;<<<<<
@read_28701_28881_324a
TGCAAGGCCGCATCGGCCAAGGCCAAGATATAGGT
+
<<<<7<<<<<<<<<<<<;6<<<:;7<<<<;<<<<<
@read_28701_28881_324b
TGCAAGGCCGCATCGGCCAAGGCCAAGATATAGGT
+
<<<<7<<<<<<<<<<<<;6<<<:;7<<<<;<<<<<
@read_28701_28881_324c
TGCAAGGCCGCATCGGCCAAGGCCAAGATATAGGT
+
<<<<7<<<<<<<<<<<<;6<<<:;7<<<<;<<<<<
@test_clipped1
AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG
+
<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<
@test_clipped1
AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG
+
<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<
"""
        )
        lane2 = mbf.align.AlignedSample(
            "test_lane2", bam_job, genome, is_paired=True, vid="AA123"
        )
        with pytest.raises(ValueError):
            lane2.to_fastq("nope.fastq")  # no support for paired end data at this point


@pytest.mark.usefixtures("new_pipegraph")
class TestQualityControl:
    def prep_lane(self):
        from mbf.sampledata import get_human_22_fake_genome

        # straight from chr22 of the human genome
        genome = get_human_22_fake_genome()

        lane = mbf.align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )
        return lane

    def _test_qc_plots(self, filename, remaining_job_count, chdir="."):
        lane = self.prep_lane()
        prune_qc(lambda job: filename in job.job_id)
        not_pruned_count = sum([1 for x in get_qc_jobs() if not x._pruned])
        assert not_pruned_count == remaining_job_count  # plot cache, plot_table, plot
        ppg.run_pipegraph()
        if chdir == "..":
            fn = lane.result_dir / chdir / filename
        else:
            fn = lane.result_dir / chdir / f"{lane.name}_{filename}"
        assert_image_equal(fn, suffix="_" + filename)

    def test_qc_complexity(self):
        self._test_qc_plots("complexity.png", 4 if is_ppg2() else 3)

    def test_qc_strandedness(self):
        self._test_qc_plots("strandedness.png", 4 if is_ppg2() else 3)

    def test_qc_reads_per_biotype(self):
        self._test_qc_plots("reads_per_biotype.png", 1)

    def test_qc_alignment_statistics(self):
        self._test_qc_plots("alignment_statistics.png", 1, "..")

    def test_qc_subchromal_distribution(self):
        self._test_qc_plots("subchromosomal_distribution.png", 4 if is_ppg2() else 3)

    def test_qc_splice_sites(self):
        self._test_qc_plots("splice_sites.png", 4 if is_ppg2() else 3)

    def test_alignment_stats(self):
        from mbf.sampledata import get_human_22_fake_genome

        genome = get_human_22_fake_genome()
        lane = mbf.align.AlignedSample(
            "test_lane",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
        )  # index creation is automatic
        counts = {"get_bam": 0}

        def get_bam():
            counts["get_bam"] += 1

            class DummySam:
                mapped = 5
                unmapped = 10

                def __enter__(self):
                    return self

                def __exit__(self, *args):
                    pass

            return DummySam()

        lane.get_bam = get_bam
        assert lane.get_alignment_stats() == {"Mapped": 5, "Unmapped": 10}
        assert counts["get_bam"] == 1

        class DummyAlignerWithout:
            pass

        lane = mbf.align.AlignedSample(
            "test_lane2",
            get_sample_data(Path("mbf_align/rnaseq_spliced_chr22.bam")),
            genome,
            False,
            "AA123",
            aligner=DummyAlignerWithout(),
        )  # index creation is automatic
        lane.get_bam = get_bam
        assert counts["get_bam"] == 1
        assert lane.get_alignment_stats() == {"Mapped": 5, "Unmapped": 10}
        assert counts["get_bam"] == 2

        class DummyAlignerWith:
            def get_alignment_stats(self, bam_filename):
                assert (
                    Path(bam_filename).resolve()
                    == get_sample_path("mbf_align/rnaseq_spliced_chr22.bam").resolve()
                )
                return {"Hello": 23}

        lane = mbf.align.AlignedSample(
            "test_lane3",
            get_sample_data("mbf_align/rnaseq_spliced_chr22.bam"),
            genome,
            False,
            "AA123",
            aligner=DummyAlignerWith(),
        )  # index creation is automatic
        lane.get_bam = get_bam
        assert counts["get_bam"] == 2
        assert lane.get_alignment_stats() == {"Hello": 23}
        assert counts["get_bam"] == 2

    def test_chromosome_mapping(self):
        bam_path = get_sample_data(Path("mbf_align/ex2.bam"))
        bam_job = ppg.FileInvariant(bam_path)
        genome = DummyGenome()
        lane = mbf.align.AlignedSample("test_lane", bam_job, genome, False, "AA123")
        assert lane.name == "test_lane"
        assert lane.load()[0] is bam_job
        assert isinstance(lane.load()[1], ppg.FileInvariant)
        assert lane.genome is genome
        assert not lane.is_paired
        assert lane.vid == "AA123"

        with pytest.raises(ValueError):
            mbf.align.AlignedSample("test_lane", bam_job, genome, False, "AA123")
        lane2 = mbf.align.AlignedSample("test_lane2", bam_job, genome, True, "AA123")
        assert lane2.is_paired

        b = lane.get_bam()
        assert isinstance(b, pysam.Samfile)
        b
