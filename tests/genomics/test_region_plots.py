import pytest
import pandas as pd
import pypipegraph2 as ppg2

from mbf.genomics.regions import GenomicRegions
from mbf.genomics.regions import plots
from mbf.align.lanes import AlignedSample
from mbf.sampledata import get_human_22_fake_genome, get_sample_path
from mbf.qualitycontrol.testing import assert_image_equal
import mbf_bam


@pytest.mark.usefixtures("new_pipegraph_no_qc")
class TestPlotAveragedCoverage:
    def test_basic(self):
        genome = get_human_22_fake_genome()

        sample = AlignedSample(
            "sampleA",
            mbf_bam.job_reheader_and_rename_chromosomes(
                get_sample_path("mbf_align/chipseq_chr22.bam"),
                "22.bam",
                {"chr22": "22"},
            ),
            genome,
            is_paired=False,
            vid=None,
        )
        intervals = []
        k = 1000
        for ii in range(int(16.5e6), int(16.5e6) + 500_000, k):
            intervals.append({"chr": "22", "start": ii, "stop": ii + k})

        gr = GenomicRegions("test", lambda: pd.DataFrame(intervals), [], genome)

        p = plots.PlotAveragedCoverage(gr, None)
        p.add_sample(sample)
        p.plot_options(title="example plot %i regions", alpha=1)
        p.render("basic.pdf")
        ppg2.run()
        assert_image_equal("basic.pdf")

    def test_substract(self):
        genome = get_human_22_fake_genome()

        sample = AlignedSample(
            "sampleA",
            mbf_bam.job_reheader_and_rename_chromosomes(
                get_sample_path("mbf_align/chipseq_chr22.bam"),
                "22.bam",
                {"chr22": "22"},
            ),
            genome,
            is_paired=False,
            vid=None,
        )
        intervals = []
        k = 1000
        for ii in range(int(16.5e6), int(16.5e6) + 500_000, k):
            intervals.append({"chr": "22", "start": ii, "stop": ii + k})

        gr = GenomicRegions("test", lambda: pd.DataFrame(intervals), [], genome)

        p = plots.PlotAveragedCoverage(gr, None)
        p.add_sample(sample, "shu", color="red", normalize=sample)
        p.plot_options(
            title="example plot %i regions",
            alpha=1,
            y_scale_args={"name": "avg signal normed to respective IgG"},
        )
        p.render("basic.pdf")
        ppg2.run()
        assert_image_equal("basic.pdf")
