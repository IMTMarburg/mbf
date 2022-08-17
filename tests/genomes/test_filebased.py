import pytest
from pathlib import Path
import pypipegraph as ppg
from mbf.genomes import FileBasedGenome, InteractiveFileBasedGenome
from mbf.genomes.common import iter_fasta, ProkaryoticCode
from mbf.externals.util import UpstreamChangedError
from pandas.testing import assert_frame_equal
from mbf.sampledata import get_sample_data


@pytest.mark.usefixtures("new_pipegraph")
class TestFilebased:
    def test_fasta_and_gtf_indexing(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.download_genome()
        g.job_transcripts()
        g.job_genes()
        with pytest.raises(ValueError):
            g.df_transcripts
        with pytest.raises(ValueError):
            g.df_genes
        ppg.run_pipegraph()
        assert g.find_file("genome.fasta").exists()
        assert g.find_prebuild("genome.fasta") == g.genome_fasta_dependencies
        assert g.find_file("genome.fasta").with_suffix(".fasta.fai").exists()
        for should_file, actual_file in [
            (
                get_sample_data(
                    "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
                ),
                g.find_file("genome.fasta"),
            ),
            (
                get_sample_data(
                    "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
                ),
                g.find_file("cdna.fasta"),
            ),
        ]:
            should = dict(iter_fasta(should_file))
            should = {k[: k.find(b" ")]: v for (k, v) in should.items()}
            actual = dict(iter_fasta(actual_file))
            if should != actual:
                assert len(should) == len(actual)
                assert set(should.keys()) == set(actual.keys())
                assert False == should_file  # noqa:E712
        tf = g.df_transcripts
        assert "BAF35033" in tf.index
        assert not hasattr(g, "_transcripts")
        assert tf.loc["BAF35033"].exons == ((1313, 2816),)

        gf = g.df_genes
        assert len(gf) == 246
        # transcript_stable_ids is tuples, this genome has only one transcript
        # per gene
        assert set([len(x) for x in gf.transcript_stable_ids]) == set([1])

        assert g.find_file("pep.fasta").exists()
        assert g.find_prebuild("pep.fasta") == g.protein_fasta_dependencies
        assert g.find_file("pep.fasta").with_suffix(".fasta.fai").exists()
        assert (
            g.get_protein_sequence("BAF35037")
            == "MFKFINRFLNLKKRYFYIFLINFFYFFNKCNFIKKKKIYKKIITKKFENYLLKLIIQKYAK"
        )

    def test_build_index(self):
        from mbf.externals.aligners.subread import Subread

        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.download_genome()
        subread = Subread()
        index = g.build_index(subread)
        ppg.run_pipegraph()
        assert len(g.df_transcripts) > 0
        assert len(g.get_gtf()) > 0
        assert len(g.df_genes) > 0
        assert len(g.df_proteins) > 0
        assert (Path(index.filenames[0]).parent / "subread_index.reads").exists()

    def test_cdna_creation(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            None,
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.download_genome()
        g.job_transcripts()
        ppg.run_pipegraph()
        assert g.find_file("genome.fasta").exists()
        assert g.find_file("genome.fasta").with_suffix(".fasta.fai").exists()
        tf = g.df_transcripts
        assert "BAF35033" in tf.index
        assert tf.loc["BAF35033"].exons == ((1313, 2816),)

        should = dict(
            iter_fasta(
                get_sample_data(
                    "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
                )
            )
        )
        should = {k[: k.find(b" ")]: v for (k, v) in should.items()}
        actual = dict(iter_fasta(g.find_file("cdna.fasta")))
        if actual != should:
            assert not set(should.keys()).difference(
                set(actual.keys())
            )  # they are all here, we just have more (tRNA...)
            for k in should:
                assert actual[k] == should[k]

    def test_empty_gtf_and_cdna_and_protein(self):
        from mbf.externals.aligners.subread import Subread

        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            None,
            None,
        )
        g.download_genome()
        assert g.gtf_filename is None
        assert g.cdna_fasta_filename is None
        g.job_transcripts()
        g.job_genes()
        g.job_proteins()
        subread = Subread()
        index = g.build_index(subread)
        ppg.run_pipegraph()
        assert len(g.df_transcripts) == 0
        assert len(g.get_gtf()) == 0
        assert len(g.df_genes) == 0
        assert len(g.df_proteins) == 0
        assert (Path(index.filenames[0]).parent / "subread_index.reads").exists()

    def test_protein_creation(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            None,
            ProkaryoticCode(),
        )
        g.download_genome()
        g.job_transcripts()
        ppg.run_pipegraph()

        should = dict(
            iter_fasta(
                get_sample_data(
                    "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
                )
            )
        )
        should = {k[: k.find(b" ")]: v for (k, v) in should.items()}
        actual = dict(iter_fasta(g.find_file("pep.fasta")))
        if actual != should:
            assert not set(should.keys()).difference(
                set(actual.keys())
            )  # they are all here, we just have more (tRNA...)
            for k in should:
                if actual[k] != should[k]:
                    print(k)
                    print(len(actual[k]))
                    print(len(should[k]))

                    print(actual[k])
                    print(should[k])
                    # print(g.get_cds_sequence(k.decode('utf-8')))
                # else:
                # print('ok', k)
                # assert actual[k] == should[k]
            assert False

    def test_job_creating_fasta(self, new_pipegraph):
        new_pipegraph.quiet = False

        def gen_fasta():
            import shutil

            shutil.copy(
                get_sample_data(
                    "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
                ),
                "shu.fasta.gz",
            )

        fasta_job = ppg.FileGeneratingJob("shu.fasta.gz", gen_fasta)
        g = FileBasedGenome(
            "Candidatus_carsonella",
            fasta_job,
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            None,
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.download_genome()
        ppg.run_pipegraph()
        assert (
            g.get_cdna_sequence("BAF35032")
            == "ATGAATACTATATTTTCAAGAATAACACCATTAGGAAATGGTACGTTATGTGTTATAAGAAT"
            "TTCTGGAAAAAATGTAAAATTTTTAATACAAAAAATTGTAAAAAAAAATATAAAAGAAAAAATAG"
            "CTACTTTTTCTAAATTATTTTTAGATAAAGAATGTGTAGATTATGCAATGATTATTTTTTTTAAA"
            "AAACCAAATACGTTCACTGGAGAAGATATAATCGAATTTCATATTCACAATAATGAAACTATTGT"
            "AAAAAAAATAATTAATTATTTATTATTAAATAAAGCAAGATTTGCAAAAGCTGGCGAATTTTTAG"
            "AAAGACGATATTTAAATGGAAAAATTTCTTTAATAGAATGCGAATTAATAAATAATAAAATTTTA"
            "TATGATAATGAAAATATGTTTCAATTAACAAAAAATTCTGAAAAAAAAATATTTTTATGTATAAT"
            "TAAAAATTTAAAATTTAAAATAAATTCTTTAATAATTTGTATTGAAATCGCAAATTTTAATTTTA"
            "GTTTTTTTTTTTTTAATGATTTTTTATTTATAAAATATACATTTAAAAAACTATTAAAACTTTTA"
            "AAAATATTAATTGATAAAATAACTGTTATAAATTATTTAAAAAAGAATTTCACAATAATGATATT"
            "AGGTAGAAGAAATGTAGGAAAGTCTACTTTATTTAATAAAATATGTGCACAATATGACTCGATTG"
            "TAACTAATATTCCTGGTACTACAAAAAATATTATATCAAAAAAAATAAAAATTTTATCTAAAAAA"
            "ATAAAAATGATGGATACAGCAGGATTAAAAATTAGAACTAAAAATTTAATTGAAAAAATTGGAAT"
            "TATTAAAAATATAAATAAAATTTATCAAGGAAATTTAATTTTGTATATGATTGATAAATTTAATA"
            "TTAAAAATATATTTTTTAACATTCCAATAGATTTTATTGATAAAATTAAATTAAATGAATTAATA"
            "ATTTTAGTTAACAAATCAGATATTTTAGGAAAAGAAGAAGGAGTTTTTAAAATAAAAAATATATT"
            "AATAATTTTAATTTCTTCTAAAAATGGAACTTTTATAAAAAATTTAAAATGTTTTATTAATAAAA"
            "TCGTTGATAATAAAGATTTTTCTAAAAATAATTATTCTGATGTTAAAATTCTATTTAATAAATTT"
            "TCTTTTTTTTATAAAGAATTTTCATGTAACTATGATTTAGTGTTATCAAAATTAATTGATTTTCA"
            "AAAAAATATATTTAAATTAACAGGAAATTTTACTAATAAAAAAATAATAAATTCTTGTTTTAGAA"
            "ATTTTTGTATTGGTAAATGA"
        )

    def test_multiple_fasta_files(self, new_pipegraph):
        import tempfile

        tf = tempfile.NamedTemporaryFile(suffix=".fasta")
        tf.write(b">Extra\nAGTC")
        tf.flush()
        g = FileBasedGenome(
            "Candidatus_carsonella",
            [
                get_sample_data(
                    "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
                ),
                tf.name,
            ],
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.download_genome()
        ppg.run_pipegraph()
        assert g.get_genome_sequence("Extra", 0, 4) == "AGTC"
        assert g.get_chromosome_lengths() == {"Extra": 4, "Chromosome": 159662}

        # test that changing the fasta leads to an explosion
        new_pipegraph.new_pipegraph()
        tf.seek(0, 0)
        tf.write(b">Extra\nAGTCA")
        tf.flush()
        g = FileBasedGenome(
            "Candidatus_carsonella",
            [
                get_sample_data(
                    "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
                ),
                tf.name,
            ],
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            None,
        )
        g.download_genome()
        assert (
            not ">Extra\nAGTCA"
            in Path(
                "cache/FileBasedGenome/Candidatus_carsonella/dna/genome.fasta"
            ).read_text()
        )

        if hasattr(ppg, "is_ppg2"):
            ppg.run_pipegraph()
            assert (
                ">Extra\nAGTCA"
                in Path(
                    "cache/FileBasedGenome/Candidatus_carsonella/dna/genome.fasta"
                ).read_text()
            )
            # for ppg2, there is no issue with redefinition - downstreams will just rerrun
        else:
            # with ppg2, we get an acception...
            with pytest.raises(UpstreamChangedError):
                ppg.run_pipegraph()

    def test_get_gtf_using_additional_gtf(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.get_additional_gene_gtfs = lambda: [
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.additional.gtf.gz"
            )
        ]
        g.download_genome()
        j = g.job_genes()
        query = get_sample_data(
            "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.additional.gtf.gz"
        )
        for x in j.prerequisites:
            if hasattr(x, "files"):
                if (query) in x.files or Path(query) in x.files:
                    break
        else:
            assert False  # wrong preqs
        ppg.run_pipegraph()
        assert "TEST1_001" in g.df_genes.index

    def test_genes_unique_check(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.get_additional_gene_gtfs = lambda: [
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            )
        ]
        g.download_genome()
        job = g.job_genes()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert "gene_stable_ids were not unique" in str(job.exception)

    def test_transcripts_unique_check(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        g.get_additional_gene_gtfs = lambda: [
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.more_transcripts.gtf.gz"
            )
        ]
        g.download_genome()
        job = g.job_transcripts()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert "transcript_stable_ids were not unique" in str(job.exception)

    def test_genes_wrong_start_stop_order(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.broken.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        job = g.job_genes()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert "start > stop" in str(job.exception)

    def test_transcript_wrong_order(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.transcript_wrong_order.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        job = g.job_transcripts()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert "start > stop" in str(job.exception)

    def test_transcript_exon_outside_transcript(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.transcript_exon_outside.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        job = g.job_transcripts()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert "Exon outside of transcript" in str(job.exception)
        assert isinstance(job.exception, ValueError)

    def test_transcript_transcript_outside_gene(self):
        g = FileBasedGenome(
            "Candidatus_carsonella",
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.transcript_outside_gene.gtf.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
            ),
            get_sample_data(
                "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
            ),
        )
        job = g.job_transcripts()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert "Transcript outside of gene" in str(job.exception)
        assert isinstance(job.exception, ValueError)

    def test_example_genome_and_interactive(self, new_pipegraph):
        from mbf.sampledata import get_Candidatus_carsonella_ruddii_pv

        g = get_Candidatus_carsonella_ruddii_pv()
        g.download_genome()
        g.job_genes()
        ppg.run_pipegraph()
        assert g.get_chromosome_lengths()
        new_pipegraph.new_pipegraph()
        Path(g.cache_dir, "lookup", "df_genes.msgpack").unlink()  # so we rerun this
        ia = InteractiveFileBasedGenome(
            "shu",
            g.find_file("genome.fasta"),
            g.find_file("cdna.fasta"),
            g.find_file("pep.fasta"),
            g.find_file("genes.gtf"),
            g.cache_dir,
        )
        assert ia.get_chromosome_lengths() == g.get_chromosome_lengths()
        print(ia.job_genes())
        ppg.run_pipegraph()
        print(ia.df_genes)
        print(g.df_genes)
        assert_frame_equal(ia.df_genes, g.df_genes)
        ppg.util.global_pipegraph = None
        Path(g.cache_dir, "lookup", "df_genes.msgpack").unlink()  # so we rerun this
        ia2 = InteractiveFileBasedGenome(
            "shu",
            g.find_file("genome.fasta"),
            g.find_file("cdna.fasta"),
            g.find_file("pep.fasta"),
            g.find_file("genes.gtf"),
            g.cache_dir,
        )
        ia.job_genes()
        assert ia2.get_chromosome_lengths() == g.get_chromosome_lengths()
        assert_frame_equal(ia2.df_genes, g.df_genes)

    def test_get_true_chromosomes(self):
        from mbf.sampledata import get_Candidatus_carsonella_ruddii_pv

        g = get_Candidatus_carsonella_ruddii_pv()
        ppg.run_pipegraph()
        assert set(g.get_chromosome_lengths()) == set(g.get_true_chromosomes())
