from mbf.externals.util import (
    to_string,
    to_bytes,
    chmod,
    lazy_method,
    lazy_lookup,
    download_http,
    download_file_and_gzip,
)

import pytest
import requests_mock
import os
from pathlib import Path
import gzip


def test_to_string():
    a = "für".encode("utf-8")
    b = "für"
    assert to_string(b) is b
    assert to_string(a) == b


def test_to_bytes():
    a = "für".encode("utf-8")
    b = "für"
    assert to_bytes(b) == a
    assert to_bytes(a) is a


def test_chmod():
    import tempfile

    tf = tempfile.NamedTemporaryFile()
    assert not os.access(tf.name, os.X_OK)
    chmod(tf.name, 0o777)
    assert os.access(tf.name, os.X_OK)


def test_lazy_method():
    class Shu:
        def __init__(self):
            self.counter = 0

        @lazy_method
        def up(self):
            self.counter += 1
            return self.counter

    x = Shu()
    assert x.up() == 1
    assert x.up() == 1


def test_lazy_lookup():
    class Shu:
        def __init__(self):
            self.counter = 0

        @lazy_lookup
        def up(self, value):
            self.counter += 1
            return self.counter * value

    x = Shu()
    assert x.up(1) == 1
    assert x.up(1) == 1
    assert x.up(2) == 4
    assert x.up(3) == 9
    assert x.up(2) == 4
    assert x.up(3) == 9


def test_lazy_lookup_on_non_method():
    global test_lazy_lookup_counter
    test_lazy_lookup_counter = 0
    @lazy_lookup
    def up(value):
        global test_lazy_lookup_counter
        test_lazy_lookup_counter += 1
        return test_lazy_lookup_counter * value

    assert up(1) == 1
    assert up(1) == 1
    assert up(2) == 4
    assert up(3) == 9
    assert up(2) == 4
    assert up(3) == 9





def test_download_404():
    with requests_mock.Mocker() as m:
        m.get("http://test.com", text="argh", status_code=404)
        with pytest.raises(ValueError):
            download_http("http://test.com", "downloaded")


def test_download_file_and_gzip(no_pipegraph):
    should = "hello world[\n"
    with requests_mock.Mocker() as m:
        m.get("http://test.com", text=should)
        with pytest.raises(ValueError):
            download_file_and_gzip("http://test.com", "test.gz.not_gz")
        download_file_and_gzip("http://test.com", "test.gz")
        assert Path("test.gz").exists()
        with gzip.GzipFile("test.gz") as op:
            actual = op.read().decode("utf-8")
        assert actual == should
