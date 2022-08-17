import pypipegraph as ppg
from mbf.externals.peak_callers import PeakZilla, MACS2


class TestPeakzilla:
    def test_quick_run(self, new_pipegraph):
        from mbf.sampledata import get_sample_path, get_human_22_fake_genome

        from mbf.align.lanes import AlignedSample
        import mbf.qualitycontrol

        new_pipegraph.quiet = False
        mbf.qualitycontrol.disable_qc()
        input_file = get_sample_path("mbf_externals/input.bam")
        background_file = get_sample_path("mbf_externals/background.bam")
        genome = get_human_22_fake_genome()
        input = AlignedSample("input", input_file, genome, is_paired=False, vid="AA000")
        background = AlignedSample(
            "background", background_file, genome, is_paired=False, vid="AA001"
        )

        a = PeakZilla()
        gr = a.call_peaks(input, background, {"-c": "1.01", "-s": "0.1"})
        gr.write()
        ppg.util.global_pipegraph.run()
        assert len(gr.df) == 37
        assert "AA000" in gr.vid
        assert "AA001" in gr.vid


class TestMacs2:
    def test_quick_run(self, new_pipegraph):
        from mbf.sampledata import get_sample_path, get_human_22_fake_genome

        from mbf.align.lanes import AlignedSample
        import mbf.qualitycontrol

        new_pipegraph.quiet = False
        mbf.qualitycontrol.disable_qc()
        input_file = get_sample_path("mbf_externals/input.bam")
        background_file = get_sample_path("mbf_externals/background.bam")
        genome = get_human_22_fake_genome()
        input = AlignedSample("input", input_file, genome, is_paired=False, vid="AA000")
        background = AlignedSample(
            "background", background_file, genome, is_paired=False, vid="AA001"
        )

        a = MACS2()
        gr = a.call_peaks(input, background, {"--nomodel": "", "--extsize": "147"})
        gr.write()
        ppg.util.global_pipegraph.run()
        assert len(gr.df) == 24
        assert "AA000" in gr.vid
        assert "AA001" in gr.vid
