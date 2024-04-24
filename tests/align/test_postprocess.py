from pathlib import Path
import pypipegraph as ppg
from mbf.align import Sample, AlignedSample
from mbf.sampledata import get_sample_path, get_human_22_fake_genome
import mbf.qualitycontrol


def test_align_and_extract_umis(new_pipegraph):
    from mbf.align.post_process import AnnotateFastqBarcodes

    for folder in [
        get_sample_path(Path("mbf_align/sample_extract_barcodes")),
        get_sample_path(Path("mbf_align/sample_extract_barcodes_gz")),
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


def test_add_chr(new_pipegraph):
    input = mbf.align.AlignedSample(
        "tesst",
        get_sample_path(Path("mbf_align/one_read_per_chr.bam")),
        get_human_22_fake_genome(),
        False,
        "AA123",
    )
    output = input.post_process(mbf.align.post_process.AddChr())
    mbf.qualitycontrol.prune_qc(lambda _: False)

    ppg.run_pipegraph()
    inp = input.get_bam()
    input_total = sum((1 for _ in inp.fetch(until_eof=True)))
    of = output.get_bam()
    out_total = sum((1 for _ in of.fetch(until_eof=True)))
    assert out_total == input_total
    assert len(inp.references) == len(of.references)
    for chr in inp.references:
        assert chr in of.references or (len(chr) < 3 and ("chr" + chr) in of.references)
        input_count = sum((1 for _ in inp.fetch(chr)))
        output_count = sum((1 for _ in of.fetch("chr" + chr if len(chr) < 3 else chr)))
        assert input_count == output_count
        if input_count > 0:
            for ip, op in zip(
                inp.fetch(chr), of.fetch("chr" + chr if len(chr) < 3 else chr)
            ):
                assert ip.pos == op.pos
                assert ip.seq == op.seq
                assert ip.qual == op.qual
                assert ip.cigarstring == op.cigarstring


def test_add_chr_and_filter(new_pipegraph):
    input = mbf.align.AlignedSample(
        "tesst",
        get_sample_path(Path("mbf_align/one_read_per_chr.bam")),
        get_human_22_fake_genome(),
        False,
        "AA123",
    )
    accepted_chrs = ["1", "KI270392.1", "nosuchchrom"]
    output = input.post_process(mbf.align.post_process.AddChrAndFilter(accepted_chrs))
    mbf.qualitycontrol.prune_qc(lambda _: False)

    ppg.run_pipegraph()
    inp = input.get_bam()
    input_total = sum((1 for _ in inp.fetch(until_eof=True)))
    of = output.get_bam()
    out_total = sum((1 for _ in of.fetch(until_eof=True)))
    assert input_total > out_total
    assert out_total == 2
    assert len(of.references) == 2
    for chr in inp.references:
        if chr in accepted_chrs:
            assert chr in of.references or (
                len(chr) < 3 and ("chr" + chr) in of.references
            )
            input_count = sum((1 for _ in inp.fetch(chr)))
            output_count = sum(
                (1 for _ in of.fetch("chr" + chr if len(chr) < 3 else chr))
            )
            assert input_count == output_count
            if input_count > 0:
                for ip, op in zip(
                    inp.fetch(chr), of.fetch("chr" + chr if len(chr) < 3 else chr)
                ):
                    assert ip.pos == op.pos
                    assert ip.seq == op.seq
                    assert ip.qual == op.qual
                    assert ip.cigarstring == op.cigarstring
