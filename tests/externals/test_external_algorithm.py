from pathlib import Path
import pypipegraph as ppg
import pytest
from mbf.externals import ExternalAlgorithm
from mbf.externals.util import Version


class WhateverAlgorithm(ExternalAlgorithm):
    @property
    def name(self):
        return "whatever"

    @property
    def primary_binary(self):
        return "bash"

    def build_cmd(self, output_directory, ncores, return_code):
        return [
            self.primary_binary,
            "-c",
            f"echo was {return_code}; exit {return_code}",
        ]

    def get_latest_version(self):
        return "0.1"


class TestExternalStore:
    def test_passing_arguments(self, new_pipegraph):
        algo = WhateverAlgorithm()
        job = algo.run(new_pipegraph.result_dir / "whatever_output", 0)
        assert job.cores_needed == 1
        ppg.util.global_pipegraph.run()
        assert Path(job.filenames[0]).exists()
        assert (Path(job.filenames[0]).parent / "stdout.txt").read_text() == "was 0\n"
        assert (Path(job.filenames[0]).parent / "stderr.txt").read_text() == ""
        assert (Path(job.filenames[0]).parent / "cmd.txt").read_text() == (
            "bash -c echo was 0; exit 0"
        )

    def test_passing_arguments_and_returncode_issues(self, new_pipegraph):
        algo = WhateverAlgorithm()
        job = algo.run(new_pipegraph.result_dir / "mhatever_output", 1)
        with pytest.raises(ppg.RuntimeError):
            ppg.util.global_pipegraph.run()
        assert not Path(job.filenames[0]).exists()
        assert (Path(job.filenames[0]).parent / "stdout.txt").read_text() == "was 1\n"
        assert (Path(job.filenames[0]).parent / "stderr.txt").read_text() == ""


class TestUtils:
    def test_get_page(self):
        from mbf.externals.util import get_page

        assert "gtf" in get_page("http://ftp.ensembl.org/pub/release-77/")
        assert "gtf" in get_page("ftp://ftp.ensembl.org/pub/release-77/")
        assert "gtf" in get_page("ftp://ftp.ensembl.org/pub/release-77/README")
        with pytest.raises(ValueError):
            get_page("ftp://ftp.ensembl.org/pub/release-77/doesnotexist")

    def test_download_file_and_gunzip(self, new_pipegraph):
        from mbf.externals.util import download_file_and_gunzip

        download_file_and_gunzip(
            "http://ftp.ensembl.org/pub/release-77/mysql/ailuropoda_melanoleuca_core_77_1/map.txt.gz",
            "map.txt",
        )
        assert Path("map.txt").read_text() == ""

    def test_download_file_with_filename_raises(self):
        from mbf.externals.util import download_file

        with pytest.raises(ValueError):
            download_file("http://ftp.ensembl.org", "out.file")

    def test_compare_versions(self):
        assert Version("") == ""
        assert Version("") == Version("")
        assert Version("0.1") == Version("0.1")
        assert Version("0.2") > "0.1"
        assert Version("0.2") < "0.3"
        assert Version("1.5.0") < "1.6.0"
        assert Version("1.4.3") < "1.5.0"
        assert Version("1.4.3-p1") < "1.5.99"
        assert not (Version("0.4") > "0.4")
        assert not (Version("0.4") < "0.4")
        assert Version("0.4") >= "0.4"
        assert Version("0.4") <= "0.4"
        assert Version("0.4") < "0.5"
        assert Version("0.4") < "0.6"
        assert Version("1.5.0") < "1.6"
        assert str(Version("1.5")) == "1.5"
        assert repr(Version("1.5")) == 'Version("1.5")'
        assert Version("1.5.0") < Version("1.6")
        assert Version("1.6.0") > Version("1.5.99.shu")
