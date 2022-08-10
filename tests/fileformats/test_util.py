from mbf.fileformats.util import chunkify, open_file


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
