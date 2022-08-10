from pathlib import Path
import pytest
import pypipegraph as ppg
from mbf.externals.aligners.subread import Subread
from mbf.externals.aligners.star import STAR
from mbf.externals.aligners.bowtie import Bowtie


class TestSubread:
    def test_build_and_align(self, new_pipegraph):
        new_pipegraph.quiet = False
        s = Subread()
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("subread_index_dir/srf")
        build_job = s.build_index_job([data_path / "genome.fasta"], None, index_name)
        align_job = s.align_job(
            data_path / "sample.fastq",
            None,
            build_job,
            "out/out.bam",
            {"input_type": "dna"},
        )
        align_job.depends_on(build_job)
        new_pipegraph.run()
        assert (Path("out") / "out.bam").exists()
        assert s.get_alignment_stats((Path("out") / "out.bam")) == {
            "Uniquely mapped": 1,
            "Multi-mapping": 0,
            "Unmapped": 0,
        }

    def test_build_index(self, new_pipegraph):
        s = Subread()
        data_path = Path(__file__).parent / "sample_data"
        s.build_index([data_path / "genome.fasta"], None, "shu")
        assert Path("shu/stdout.txt").exists()
        assert Path("shu/subread_index.reads").exists()

    def test_raises_on_invalid_input_type(self, new_pipegraph):
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("subread_index_dir/srf")
        s = Subread()
        index_job = s.build_index_job([data_path / "genome.fasta"], None, index_name)
        with pytest.raises(ValueError):
            s.align_job(data_path / "sample.fastq", None, index_job, "out/out.bam", {})
        with pytest.raises(ValueError):
            s.align_job(
                data_path / "sample.fastq",
                None,
                index_job,
                "out/out.bam",
                {"input_type": "shu"},
            )

    def test_cant_call_subread_run(self, new_pipegraph):
        s = Subread()
        s.run("out", None)
        with pytest.raises(ppg.RuntimeError):
            new_pipegraph.run()

    def test_cant_call_subread_run2(self, new_pipegraph):
        s = Subread()
        s.run("out", ["subread-align", "something"])
        with pytest.raises(ppg.RuntimeError):
            new_pipegraph.run()

    def test_build_and_align_paired_end(self, new_pipegraph):
        import pypipegraph2 as ppg2

        new_pipegraph.quiet = False
        s = Subread()
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("subread_index_dir/srf")
        build_job = s.build_index_job(data_path / "genome.fasta", None, index_name)
        new_pipegraph.run()

        print(id(ppg2.global_pipegraph))
        new_pipegraph.new_pipegraph()  # so we can test 'pass index as path'
        print(id(ppg2.global_pipegraph))
        align_job = s.align_job(
            data_path / "sample_R1_.fastq",
            data_path / "sample_R2_.fastq",
            index_name,
            "out/out.bam",
            {"input_type": "rna"},
        )
        ppg.util.global_pipegraph.run()
        print(align_job)
        assert (Path("out") / "out.bam").exists()

    def test_get_index_version_range(self, new_pipegraph):
        s = Subread()
        print(s.version)
        assert s.version.count(".") == 2  # heuristic of course
        s.version = "1.4.3-p1"
        assert s.get_index_version_range() == ("0.1", "1.5.99")
        s.version = "1.6.3"
        assert s.get_index_version_range() == ("1.6", None)


class TestSTAR:
    def test_build_and_align(self, new_pipegraph):
        new_pipegraph.quiet = False
        s = STAR()
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("star/srf")
        build_job = s.build_index_job(
            data_path / "genome.fasta", data_path / "genes.gtf", index_name
        )
        align_job = s.align_job(
            data_path / "sample.fastq",
            None,
            build_job,
            "out/out.bam",
            parameters={"--runRNGseed": "5555"},
        )
        new_pipegraph.run()
        assert (Path("out") / "out.bam").exists()
        assert "--runRNGseed 5555" in (Path("out") / "cmd.txt").read_text()
        assert s.get_alignment_stats((Path("out") / "out.bam")) == {
            "Number of reads mapped to multiple loci": 0,
            "Number of reads mapped to too many loci": 0,
            "Uniquely mapped reads number": 1,
            "Unmapped": 0,
        }

    def test_build_and_align_paired_end(self, new_pipegraph):
        new_pipegraph.quiet = False
        s = STAR()
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("star/srf")
        build_job = s.build_index_job(
            [data_path / "genome.fasta"], data_path / "genes.gtf", index_name
        )
        align_job = s.align_job(
            data_path / "sample_R1_.fastq",
            data_path / "sample_R2_.fastq",
            build_job,
            "out/out.bam",
            parameters={"--runRNGseed": "5555"},
        )
        align_job.depends_on(build_job)
        new_pipegraph.run()
        assert (Path("out") / "out.bam").exists()
        assert "--runRNGseed 5555" in (Path("out") / "cmd.txt").read_text()

    def test_cant_call_star_run(self, new_pipegraph):
        s = STAR()
        s.run("out", None)
        with pytest.raises(ppg.RuntimeError):
            new_pipegraph.run()

    def test_build_raises_on_multiple_fasta(self, new_pipegraph):
        s = STAR()
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("star/srf")

        with pytest.raises(ValueError):
            s.build_index_job(
                [data_path / "genome.fasta", data_path / "genome.fasta"],
                data_path / "genes.gtf",
                index_name,
            )
        with pytest.raises(ValueError):
            s.build_index_job(data_path / "genome.fasta", None, index_name)


class TestBowtie:
    def test_build_and_align(self, new_pipegraph):
        new_pipegraph.quiet = False
        s = Bowtie()
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("bowtie/srf")
        build_job = s.build_index_job(
            data_path / "genome.fasta", data_path / "genes.gtf", index_name
        )
        align_job = s.align_job(
            data_path / "sample.fastq",
            None,
            build_job,
            "out/out.bam",
            parameters={"-k": "2"},
        )
        align_job.depends_on(build_job)
        new_pipegraph.run()
        assert (Path("out") / "out.bam").exists()
        assert "-k 2" in (Path("out") / "cmd.txt").read_text()

    def test_build_and_align_paired_end(self, new_pipegraph):
        new_pipegraph.quiet = False
        s = Bowtie()
        data_path = Path(__file__).parent / "sample_data"
        index_name = Path("bowtie/srf")
        build_job = s.build_index_job(
            [data_path / "genome.fasta"], data_path / "genes.gtf", index_name
        )
        align_job = s.align_job(
            data_path / "sample_R1_.fastq",
            data_path / "sample_R2_.fastq",
            build_job,
            "out/out.bam",
            parameters={},
        )
        align_job.depends_on(build_job)
        new_pipegraph.run()
        assert (Path("out") / "out.bam").exists()
