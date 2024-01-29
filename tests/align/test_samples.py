import pytest
import gzip
from pathlib import Path
from pypipegraph import FileGeneratingJob, MultiFileGeneratingJob
import requests_mock
import pypipegraph as ppg

from mbf.align import (
    FASTQsFromFile,
    FASTQsFromFiles,
    FASTQsFromFolder,
    FASTQsFromJob,
    FASTQsFromURLs,
    FASTQsFromAccession,
    FASTQsFromPrefix,
    build_fastq_strategy,
    FASTQsFromMRNAs,
    FASTQsJoin,
)
from mbf.align import Sample, PairingError, fastq2, sanity_check
from mbf.align._common import read_fastq_iterator
from mbf.sampledata import get_sample_data
import attr


def test_FASTQsFromFile():
    fn = Path(
        get_sample_data(Path("mbf_align/sample_a") / ".." / "sample_a" / "a.fastq")
    )
    o = FASTQsFromFile(fn)
    assert o() == [(fn.resolve(),)]


def test_FASTQsFromFileRaisesOnMissing():
    fn = get_sample_data(Path("mbf_align/sample_a") / "a.fastq.nosuchfile")
    with pytest.raises(IOError):
        FASTQsFromFile(fn)


def test_FASTQsFromFilePaired():
    fn = Path(get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"))
    fn2 = Path(
        get_sample_data(
            Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R2_.fastq.gz"
        )
    )
    o = FASTQsFromFile(fn, fn2)
    assert o() == [(fn.resolve(), fn2.resolve())]


def test_FASTQsFromFilePairedMissingR2():
    fn = get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz")
    fn2 = get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz.nosuchfile")
    with pytest.raises(IOError):
        FASTQsFromFile(fn, fn2)


def test_FASTQsFromFilesPaired():
    fn = Path(get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"))
    fn2 = Path(
        get_sample_data(
            Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R2_.fastq.gz"
        )
    )
    o = FASTQsFromFiles([fn, fn2])
    assert o() == [(fn.resolve(), fn2.resolve())]


def test_FASTQsFromFilesPaired_build_strategy():
    fn = Path(get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"))
    fn2 = Path(
        get_sample_data(
            Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R2_.fastq.gz"
        )
    )
    o = build_fastq_strategy([fn, fn2])
    assert o() == [(fn.resolve(), fn2.resolve())]


def test_FASTQsFromFolder():
    folder = Path(get_sample_data(Path("mbf_align/sample_a")))
    o = FASTQsFromFolder(folder)
    import pprint

    pprint.pprint(o())
    assert o() == [
        ((folder / "a.fastq").resolve(),),
        ((folder / "b.fastq.gz").resolve(),),
    ]


def test_fastqs_join():
    fn = Path(get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"))
    fn2 = Path(
        get_sample_data(
            Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R2_.fastq.gz"
        )
    )
    a = FASTQsFromFiles([fn, fn2])
    b = FASTQsFromFile(fn)
    c = FASTQsFromFile(fn2)
    d = FASTQsJoin([a, b, c])
    o = d()
    assert o == [(fn.resolve(), fn2.resolve()), (fn.resolve(),), (fn2.resolve(),)]


def test_FASTQsFromFolder_raises_on_non_existing():
    with pytest.raises(IOError):
        FASTQsFromFolder("shu")


def test_FASTQsFromFolder_raises_on_no_fastqs():
    with pytest.raises(ValueError):
        FASTQsFromFolder(get_sample_data(Path("mbf_align/sample_f")))


def test_FASTQsFromFolder_raises_on_not_a_folder():
    with pytest.raises(ValueError):
        FASTQsFromFolder(get_sample_data(Path("mbf_align/sample_a") / "a.fastq"))


def test_FASTQsFromFolderPaired():
    folder = Path(get_sample_data(Path("mbf_align/sample_b")))
    o = FASTQsFromFolder(folder)
    assert o() == [
        ((folder / "a_R1_.fastq.gz").resolve(), (folder / "a_R2_.fastq.gz").resolve())
    ]


def test_FASTQsFromFolderR2_but_missing_any_r1():
    folder = get_sample_data(Path("mbf_align/sample_c"))
    o = FASTQsFromFolder(folder)
    with pytest.raises(ValueError):
        o()


def test_FASTQsFromFolder_pairing_files_fails():
    folder = get_sample_data(Path("mbf_align/sample_d"))
    o = FASTQsFromFolder(folder)
    with pytest.raises(ValueError):
        o()


def test_FASTQsFromFolder_pairing_files_fails2():
    folder = get_sample_data(Path("mbf_align/sample_e"))
    o = FASTQsFromFolder(folder)
    with pytest.raises(ValueError):
        o()


def test_FASTQsFromPrefix():
    fn1 = Path(get_sample_data(Path("mbf_align/sample_d") / "a_R1_.fastq.gz"))
    fn2 = Path(get_sample_data(Path("mbf_align/sample_d") / "a_R2_.fastq.gz"))
    fn_prefix = Path(get_sample_data(Path("mbf_align/sample_d") / "a"))
    o = FASTQsFromPrefix(fn_prefix)
    str(o)
    assert o() == [(fn1.resolve(), fn2.resolve())]


def test_FASTQsFromPrefix_raises_on_non_existant():
    fn_prefix = Path("/shu/sha")
    with pytest.raises(IOError):
        FASTQsFromPrefix(fn_prefix)


def test_FASTQsFromPrefix_raises_on_non_found():
    fn_prefix = Path(get_sample_data(Path("mbf_align/sample_d") / "x"))
    with pytest.raises(ValueError):
        FASTQsFromPrefix(fn_prefix)


@pytest.mark.usefixtures("new_pipegraph_no_qc")
class TestSamples:
    def test_FASTQsFromJob(self):
        job = FileGeneratingJob("test.fastq.gz", lambda of: None)
        o = FASTQsFromJob(job)
        assert o() == [(Path("test.fastq.gz").resolve(),)]

    def test_FASTQsFromJob_R1_ok(self):
        job = FileGeneratingJob("test_R1_.fastq.gz", lambda of: None)
        o = FASTQsFromJob(job)
        assert o() == [(Path("test_R1_.fastq.gz").resolve(),)]

    def test_FASTQsFromJob_Multiple_Unpaired(self):
        job = MultiFileGeneratingJob(
            ["test.fastq.gz", "test2.fastq.gz"], lambda of: None
        )
        o = FASTQsFromJob(job)
        assert o() == [
            (Path("test.fastq.gz").resolve(),),
            (Path("test2.fastq.gz").resolve(),),
        ]

    def test_FASTQsFromJob_Multiple_Unpaired_R1(self):
        job = MultiFileGeneratingJob(
            ["test_R1_.fastq.gz", "test2_R1_.fastq.gz"], lambda of: None
        )
        o = FASTQsFromJob(job)
        assert o() == [
            # 2 sorts before _
            (Path("test2_R1_.fastq.gz").resolve(),),
            (Path("test_R1_.fastq.gz").resolve(),),
        ]

    def test_FASTQsFromJob_Multiple_Paired(self):
        job = MultiFileGeneratingJob(
            [
                "test_R1_.fastq.gz",
                "test2_R1_.fastq.gz",
                "test_R2_.fastq.gz",
                "test2_R2_.fastq.gz",
            ],
            lambda of: None,
        )
        o = FASTQsFromJob(job)
        assert set(o()) == set(
            [
                # 2 sorts before _
                (
                    Path("test2_R1_.fastq.gz").resolve(),
                    Path("test2_R2_.fastq.gz").resolve(),
                ),
                (
                    Path("test_R1_.fastq.gz").resolve(),
                    Path("test_R2_.fastq.gz").resolve(),
                ),
            ]
        )

    def test_build_fastq_strategy(self):
        # single filename
        assert build_fastq_strategy(
            get_sample_data(
                Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R1_.fastq.gz"
            )
        )() == [
            (
                (
                    Path(get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"))
                ).resolve(),
            )
        ]
        assert build_fastq_strategy(
            str(
                get_sample_data(
                    Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R1_.fastq.gz"
                )
            )
        )() == [
            (
                (
                    Path(get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"))
                ).resolve(),
            )
        ]
        # multiple files - end up being paired!
        assert build_fastq_strategy(
            [
                get_sample_data(
                    Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R1_.fastq.gz"
                ),
                get_sample_data(
                    Path("mbf_align/sample_b") / ".." / "sample_b" / "a_R2_.fastq.gz"
                ),
            ]
        )() == [
            (
                Path(
                    get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz")
                ).resolve(),
                Path(
                    get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz")
                ).resolve(),
            )
        ]
        # folder
        assert build_fastq_strategy(
            get_sample_data(Path("mbf_align/sample_b") / ".." / "sample_b")
        )() == [
            (
                Path(
                    get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz")
                ).resolve(),
                Path(
                    get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz")
                ).resolve(),
            )
        ]
        # job
        assert isinstance(
            build_fastq_strategy(FileGeneratingJob("test.fastq", lambda of: None)),
            FASTQsFromJob,
        )

        # pass through
        fn = get_sample_data(Path("mbf_align/sample_a") / ".." / "sample_a" / "a.fastq")
        o = FASTQsFromFile(fn)
        assert build_fastq_strategy(o) is o

        with pytest.raises(ValueError):
            build_fastq_strategy(55)

    def test_lane(self):
        lane = Sample(
            "Sample_a", get_sample_data(Path("mbf_align/sample_a")), False, vid="VA000"
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert Path(real_job.filenames[0]).exists()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            lines = op.readlines()
            assert len(lines) == 20 + 20

    def test_paired_modes(self):
        with pytest.raises(PairingError):
            lane = Sample(
                "Sample_a",
                get_sample_data(Path("mbf_align/sample_b")),
                False,
                vid="VA000",
            )
            lane.prepare_input()

    def test_lane_paired_straight(self):
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_b")),
            False,
            vid="VA000",
            pairing="paired",
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert not Path(temp_job.filenames[1]).exists()
        assert Path(real_job.filenames[0]).exists()
        assert Path(real_job.filenames[1]).exists()
        assert "_R1_" in str(real_job.filenames[0])
        assert "_R2_" in str(real_job.filenames[1])
        assert ".fastq.gz" in str(real_job.filenames[0])
        assert ".fastq.gz" in str(real_job.filenames[1])

        for input_fn, output_fn in zip(
            [
                (get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz")),
                (get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz")),
            ],
            real_job.filenames,
        ):
            with gzip.GzipFile(output_fn, "r") as op:
                actual = op.read()
            with gzip.GzipFile(input_fn, "r") as op:
                should = op.read()
            assert actual == should

    def test_lane_paired_filtered(self):
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_b")),
            False,
            vid="VA000",
            pairing="paired",
            fastq_processor=fastq2.Paired_Filtered(lambda *args: True),
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert not Path(temp_job.filenames[1]).exists()
        assert Path(real_job.filenames[0]).exists()
        assert Path(real_job.filenames[1]).exists()
        assert "_R1_" in str(real_job.filenames[0])
        assert "_R2_" in str(real_job.filenames[1])
        assert ".fastq.gz" in str(real_job.filenames[0])
        assert ".fastq.gz" in str(real_job.filenames[1])

        for input_fn, output_fn in zip(
            [
                (get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz")),
                (get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz")),
            ],
            real_job.filenames,
        ):
            with gzip.GzipFile(output_fn, "r") as op:
                actual = op.read()
            with gzip.GzipFile(input_fn, "r") as op:
                should = op.read()
            assert actual == should

    def test_lane_paired_as_single(self):
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_b")),
            False,
            vid="VA000",
            pairing="paired_as_single",
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert len(temp_job.filenames) == 1
        assert Path(real_job.filenames[0]).exists()
        assert len(real_job.filenames) == 1
        assert not "_R1_" in str(real_job.filenames[0])
        assert ".fastq.gz" in str(real_job.filenames[0])

        should = b""
        for input_fn in [
            (get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz")),
            (get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz")),
        ]:
            with gzip.GzipFile(input_fn, "r") as op:
                should += op.read()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            actual = op.read()
        assert actual == should

    def test_lane_paired_missing_R2(self):
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_a")),
            False,
            vid="VA000",
            pairing="paired",
        )
        with pytest.raises(PairingError):
            lane.prepare_input()

    def test_lane_paired_only_first(self):
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_b")),
            False,
            vid="VA000",
            pairing="only_first",
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert len(temp_job.filenames) == 1
        assert Path(real_job.filenames[0]).exists()
        assert len(real_job.filenames) == 1
        assert not "_R1_" in str(real_job.filenames[0])
        assert ".fastq.gz" in str(real_job.filenames[0])

        should = b""
        for input_fn in [
            (get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"))
        ]:
            with gzip.GzipFile(input_fn, "r") as op:
                should += op.read()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            actual = op.read()
        assert actual == should

    def test_lane_paired_only_second(self):
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_b")),
            False,
            vid="VA000",
            pairing="only_second",
        )
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert len(temp_job.filenames) == 1
        assert Path(real_job.filenames[0]).exists()
        assert len(real_job.filenames) == 1
        assert not "_R1_" in str(real_job.filenames[0])
        assert ".fastq.gz" in str(real_job.filenames[0])

        should = b""
        for input_fn in [
            (get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz"))
        ]:
            with gzip.GzipFile(input_fn, "r") as op:
                should += op.read()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            actual = op.read()
        assert actual == should

    def test_pairing_invalid_value(self):
        with pytest.raises(ValueError):
            Sample(
                "Sample_a",
                get_sample_data(Path("mbf_align/sample_a")),
                False,
                pairing="do_what_you_want",
            )
        with pytest.raises(ValueError):
            Sample(
                "Sample_a",
                get_sample_data(Path("mbf_align/sample_a")),
                False,
                pairing=False,
            )
        with pytest.raises(ValueError):
            Sample(
                "Sample_a",
                get_sample_data(Path("mbf_align/sample_a")),
                False,
                pairing=None,
            )
        with pytest.raises(ValueError):
            Sample(
                "Sample_a",
                get_sample_data(Path("mbf_align/sample_a")),
                False,
                pairing=[5],
            )

    def test_lane_raises_on_pe_as_se(self):
        lane = Sample("Sample_a", get_sample_data(Path("mbf_align/sample_b")), False)
        with pytest.raises(PairingError):
            lane.prepare_input()

    def test_lane_with_job_generating_fastq(self):
        def gen_fastq(fn):
            with open(fn, "wb") as op:
                op.write(b"@shu\nAGTC\n+\nasdf")

        job = FileGeneratingJob("input.fastq", gen_fastq)

        lane = Sample("Sample_a", job, False, vid="VA000")
        assert lane.vid == "VA000"
        temp_job = lane.prepare_input()
        assert job in temp_job.prerequisites
        real_job = lane.save_input()
        ppg.run_pipegraph()
        assert not Path(temp_job.filenames[0]).exists()
        assert Path(real_job.filenames[0]).exists()
        with gzip.GzipFile(real_job.filenames[0], "r") as op:
            lines = op.readlines()
        assert len(lines) == 4

    def test_align(self):
        import json
        import gzip

        class FakeGenome:
            name = "FakeGenome"

            def download_genome(self):
                return []

            def job_genes(self):
                return []

            def job_transcripts(self):
                return []

            def build_index(self, aligner, fasta_to_use=None, gtf_to_use=None):
                job = ppg.FileGeneratingJob(
                    "fake_index", lambda: Path("fake_index").write_text("hello")
                )
                job.output_path = "fake_index"
                return job

        class FakeAligner:
            name = "FakeAligner"
            version = "0.1"

            def align_job(
                self,
                input_fastq,
                paired_end_filename,
                index_basename,
                output_bam_filename,
                parameters,
            ):
                def align():
                    with open(output_bam_filename, "w") as op:
                        json.dump(
                            [
                                open(input_fastq).read(200),
                                open(paired_end_filename).read(200)
                                if paired_end_filename
                                else "",
                                str(index_basename.files[0]),
                                str(parameters),
                            ],
                            op,
                        )
                    with open(str(output_bam_filename) + ".bai", "w") as op:
                        op.write("Done")

                job = ppg.MultiFileGeneratingJob(
                    [output_bam_filename, str(output_bam_filename) + ".bai"], align
                )
                job.depends_on_params("")
                return job

        aligner = FakeAligner()
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_b")),
            False,
            vid="VA000",
            pairing="paired",
        )
        genome = FakeGenome()
        params = {"shu": 123}
        aligned_lane = lane.align(aligner, genome, params)
        ppg.run_pipegraph()
        assert Path("fake_index").exists()
        assert Path("fake_index").read_text() == "hello"
        assert str(aligned_lane.load()[0].filenames[0]).endswith(lane.name + ".bam")
        assert str(aligned_lane.load()[0].filenames[1]).endswith(lane.name + ".bam.bai")
        assert Path(aligned_lane.load()[0].filenames[0]).exists()
        with open(aligned_lane.load()[0].filenames[0]) as op:
            actual = json.load(op)
        with gzip.GzipFile(
            get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz")
        ) as op:
            should_0 = op.read(200).decode("utf-8")
        with gzip.GzipFile(
            get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz")
        ) as op:
            should_1 = op.read(200).decode("utf-8")

        assert actual[0] == should_0
        assert actual[1] == should_1
        assert actual[2] == "fake_index"
        assert actual[3] == str(params)

    def test_align_parameterDependencyChecking(self):
        class FakeGenome:
            name = "FakeGenome"

            def build_index(self, aligner, fasta_to_use=None, gtf_to_use=None):
                job = ppg.FileGeneratingJob(
                    "fake_index", lambda: Path("fake_index").write_text("hello")
                )
                job.output_path = "fake_index"
                return job

        class FakeAligner:
            name = "FakeAligner"
            version = "0.1"

            def align_job(
                self,
                input_fastq,
                paired_end_filename,
                index_basename,
                output_bam_filename,
                parameters,
            ):
                job = ppg.MultiFileGeneratingJob(
                    [output_bam_filename, str(output_bam_filename) + ".bai"], lambda: 5
                )
                # job.depends_on_params("") # that's the line we check
                return job

        aligner = FakeAligner()
        lane = Sample(
            "Sample_a",
            get_sample_data(Path("mbf_align/sample_b")),
            False,
            vid="VA000",
            pairing="paired",
        )
        genome = FakeGenome()
        params = {"shu": 123}
        with pytest.raises(ppg.JobContractError):
            lane.align(aligner, genome, params)

    def test_from_url(self):
        import requests_mock

        with requests_mock.mock() as m:
            url = "https://www.imt.uni-marburg.de/sample.fastq.gz"
            m.get(url, text="hello_world")
            o = FASTQsFromURLs(url)
            ppg.run_pipegraph()
            assert len(o.target_files) == 1
            assert len(o.dependencies) == 2
            assert Path(o.dependencies[0].filenames[0]).read_text() == "hello_world"
            assert Path(o.dependencies[0].filenames[1]).read_text() == url
            assert o() == [(Path(o.dependencies[0].filenames[0]).absolute(),)]

    def test_from_url_paired(self):
        import requests_mock

        with requests_mock.mock() as m:
            url1 = "https://www.imt.uni-marburg.de/sample_R1_.fastq.gz"
            url2 = "https://www.imt.uni-marburg.de/sample_R2_.fastq.gz"
            m.get(url1, text="hello_world1")
            m.get(url2, text="hello_world2")
            o = FASTQsFromURLs([url1, url2])
            ppg.run_pipegraph()
            assert len(o.target_files) == 2
            assert len(o.dependencies) == 3
            assert Path(o.dependencies[0].filenames[0]).read_text() == "hello_world1"
            assert "_R1_.fastq.gz" in str(o.dependencies[0].filenames[0])
            assert Path(o.dependencies[0].filenames[1]).read_text() == url1
            assert Path(o.dependencies[1].filenames[0]).read_text() == "hello_world2"
            assert Path(o.dependencies[1].filenames[1]).read_text() == url2
            assert o() == [
                (
                    Path(o.dependencies[0].filenames[0]).absolute(),
                    Path(o.dependencies[1].filenames[0]).absolute(),
                )
            ]

    def test_from_url_detects_404(self):
        with requests_mock.mock() as m:
            url = "https://www.imt.uni-marburg.de/sample.fastq.gz"
            m.get(url, text="hello_world", status_code=404)
            o = FASTQsFromURLs(url)
            o.download_files()
            with pytest.raises(ppg.RuntimeError):
                ppg.run_pipegraph()
            assert "Error return" in str(o.dependencies[0].exception)
            assert url in str(o.dependencies[0].exception)

    def test_fastqs_from_err(self):
        with requests_mock.mock() as m:
            m.get(
                "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERR2223563&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes",
                text="""run_accession	fastq_ftp	fastq_md5	fastq_bytes
ERR2223563	ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_2.fastq.gz	0e29c053bcd31072c8bed9eddece1cec;5d848b65379c195fe158a5d7324b4a18	1170312089;1246298835""",
            )
            o = FASTQsFromAccession("ERR2223563")
            print(o.urls)
            assert (
                o.urls[0]
                == "http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_1.fastq.gz"
            )
            assert (
                o.urls[1]
                == "http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR222/003/ERR2223563/ERR2223563_2.fastq.gz"
            )

    def test_fastqs_from_mrnas(self):
        @attr.s
        class Transcript:
            transcript_stable_id = attr.ib()
            mrna = attr.ib()

        class FakeGenome:
            name = "FakeGenome"

            def download_genome(self):
                return []

            def job_genes(self):
                return []

            def job_transcripts(self):
                return []

            transcripts = {
                "tr1": Transcript("gene1", "AGTC"),
                "tr2": Transcript("gene1", "ACCA"),
            }

        o = FASTQsFromMRNAs(["tr1", "tr2"], FakeGenome(), 2)
        target = o()[0][0]
        ppg.run_pipegraph()
        assert Path(target).exists()
        with open(target) as r:
            seqs, names, quals = zip(*read_fastq_iterator(r))
        assert list(seqs) == ["AG", "GT", "TC", "AC", "CC", "CA"]


@pytest.mark.usefixtures("new_pipegraph")
class TestSamplesQC:
    def test_fastqc(self):
        from mbf.qualitycontrol import get_qc_jobs

        lane = Sample(
            "Sample_a", get_sample_data(Path("mbf_align/sample_a")), False, vid="VA000"
        )
        qc_jobs = list(get_qc_jobs())
        assert len(qc_jobs) == 1
        assert "results/lanes/Sample_a/FASTQC/sentinel.txt" in str(qc_jobs[0].filenames)
        assert lane.prepare_input() in qc_jobs[0].prerequisites


@pytest.mark.usefixtures("new_pipegraph")
def test_sanity_check():
    p = Path("fastqs")
    p.mkdir()
    (p / "A_r1_1.fastq").touch()
    (p / "A_r1_2.fastq").touch()
    (p / "B_r1_1.fastq").touch()
    (p / "B_r1_2.fastq").touch()
    sampleA = Sample("A", FASTQsFromPrefix("fastqs/A"), reverse_reads=False)
    sampleB = Sample("B", FASTQsFromPrefix("fastqs/B"), reverse_reads=False)
    sanity_check({"a": sampleA, "b": sampleB})
    (p / "B_r1_3.fastq").touch()
    with pytest.raises(ValueError):
        sanity_check({"a": sampleA, "b": sampleB})
