from pathlib import Path
import pypipegraph as ppg
from mbf.align import Sample, AlignedSample
from mbf.sampledata import get_sample_path, get_human_22_fake_genome
import mbf.qualitycontrol


def test_align_and_extract_umis(new_pipegraph):
    from mbf.align.post_process import AnnotateFastqBarcodes

    for folder in [
        get_sample_path(Path("mbf.align/sample_extract_barcodes")),
        get_sample_path(Path("mbf.align/sample_extract_barcodes_gz")),
    ]:
        new_pipegraph.new_pipegraph()
        genome = get_human_22_fake_genome()

        mbf.qualitycontrol.prune_qc(lambda _: False)
        r = Sample("test", str(folder), False, pairing="only_second", vid="AA123")
        al = AlignedSample("test", str(folder / "test.bam"), genome, False, "AA123")

        x = al.post_process(AnnotateFastqBarcodes(r, {"XC": [0, 4], "XM": [7, 7 + 4]}))
        ppg.run_pipegraph()
        f = x.get_bam()
        r = next(f.fetch())
        print(r.tags)
        assert r.get_tag("XC") == "AGTC"
        assert r.get_tag("XM") == "TGAC"
