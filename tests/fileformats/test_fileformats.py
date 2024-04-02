from mbf.fileformats.util import chunkify, open_file
import pytest


def test_open_file():
    import gzip
    import tempfile
    import bz2

    tf = tempfile.TemporaryFile()
    assert open_file(tf) is tf

    tf2 = tempfile.NamedTemporaryFile(suffix=".gz", mode="w")
    g = gzip.GzipFile(tf2.name, "w")
    g.write(b"hello")
    g.close()
    tf2.flush()
    assert open_file(tf2.name).read() == b"hello"

    tf3 = tempfile.NamedTemporaryFile(suffix=".bz2", mode="w")
    b = bz2.BZ2File(tf3.name, "w")
    b.write(b"world")
    b.close()
    tf3.flush()
    assert open_file(tf3.name).read() == b"world"


def test_chunkify():
    import tempfile

    tf = tempfile.TemporaryFile("w+")
    tf.write("hello world")
    tf.flush()
    tf.seek(0, 0)
    c = list(chunkify(tf, " ", 2))
    assert c == ["hello", "world"]


def test_multi_track_bed():
    from mbf.fileformats.bed import read_bed
    import mbf.sampledata   

    fn = mbf.sampledata.get_sample_path("mbf_fileformats/multi_track.bed")
    bed = read_bed(fn)
    assert len(bed) == 3
    bed_1 = read_bed(fn, filter_to_track="input.bam_Stitched")
    assert len(bed_1) == 1
    assert bed_1[0].refseq == 'chr1'
    assert bed_1[0].position == 54272915
    assert bed_1[0].length == 54357431 - 54272915
    bed_2 = read_bed(fn, filter_to_track="Super_input.bam_Stitched")
    assert len(bed_2) == 2
    assert bed_2[0].refseq == 'chr11'
    assert bed_2[1].refseq == 'chr2'
    with pytest.raises(ValueError):
        read_bed(fn, filter_to_track="No such track")


