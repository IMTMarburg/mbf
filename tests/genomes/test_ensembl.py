import pytest
from mbf.genomes import EnsemblGenome
from mbf.externals import PrebuildManager
from mbf.externals.aligners.subread import Subread
import pypipegraph as ppg
from pypipegraph.util import checksum_file
from unittest.mock import patch


@pytest.mark.usefixtures("new_pipegraph")
class TestEnsembl:
    def test_download(self, new_pipegraph, mock_download, shared_prebuild):
        species = "Ashbya_gossypii"  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome(species, "41", prebuild_manager=shared_prebuild)

        def shorten_genome_fasta(output_path):
            with open(g.find_file("genome.fasta")) as op:
                head = op.read(1024 * 100)
            (output_path / "test.fasta").write_text(head)

        test_fasta_job = g.prebuild_manager.prebuild(
            f"ensembl/{g.species}_{g.revision}/test_fasta",
            "1",
            [],
            ["test.fasta"],
            shorten_genome_fasta,
        )
        test_fasta_job.depends_on(g.download_genome())
        g._prebuilds.append(test_fasta_job)

        subread = Subread()
        index = g.build_index(subread, "test.fasta", postfix="test")

        new_pipegraph.run()
        # note that these are not the checksums from CHECKSUMS files (those are fore
        # the gziped variants, we keep them ungziped and let the filesystem handle
        # the gzip, since we can't rely on the downstream reading gzip...
        assert (
            checksum_file(g.find_file("genome.fasta"))
            == "584a734589964a654c7c1dc23b0167ab"
        )
        assert (
            checksum_file(g.find_file("cdna.fasta"))
            == "3fc1f19ab829573169cb2488abe39211"
        )
        assert (
            checksum_file(g.find_file("genes.gtf"))
            == "8bdeec9b3db5278668dbff8b34e9d93b"
        )
        assert (
            checksum_file(g.find_file("genes.gtf"))
            == "8bdeec9b3db5278668dbff8b34e9d93b"
        )
        assert (
            checksum_file(g.find_file("pep.fasta"))
            == "9580fd44832d419c38469d657f6e2484"
        )
        with pytest.raises(OSError):
            g.find_file("no such file")
        assert index["subread_index.reads"].exists()
        assert index["subread_index.files"].exists()
        assert index["subread_index.00.b.array"].exists()

        assert g.find_file("genome.fasta.fai").exists()
        assert g.find_file("cdna.fasta.fai").exists()

        new_pipegraph.new_pipegraph()
        pb = PrebuildManager(shared_prebuild.prebuilt_path)
        g = EnsemblGenome(species, "41", prebuild_manager=pb)
        test_fasta_job = g.prebuild_manager.prebuild(
            f"ensembl/{g.species}_{g.revision}/test_fasta",
            "1",
            [],
            ["test.fasta"],
            shorten_genome_fasta,
        )
        g._prebuilds.append(test_fasta_job)
        if hasattr(ppg, "is_ppg2"):
            test_fasta_job.depends_on(g.find_prebuild("genome.fasta"))

        subread_intermediate = Subread()
        _index_intermediate = g.build_index(  # noqa: F841
            subread_intermediate, fasta_to_use="test.fasta", postfix="test"
        )
        index_genome = g.build_index(subread_intermediate)
        assert "/genome/" in str(index_genome.files[0])

        ppg.run_pipegraph()
        assert g.get_chromosome_lengths() == {
            "IV": 1_467_287,
            "MT": 23564,
            "V": 1_519_140,
            "III": 907_494,
            "II": 870_771,
            "VII": 1_800_949,
            "I": 693_414,
            "VI": 1_836_693,
        }

        assert g.get_genome_sequence("VI", 20, 30) == "ACCGCTGAGA"
        assert (
            g.get_cdna_sequence("EFAGOT00000000349")
            == "GCTCGCGTGGCGTAATGGCAACGCGTCTGACTTCTAATCAGAAGATTGTGGGTTCGACCC"
            "CCACCGTGAGTG"
        )
        assert (
            g.get_protein_sequence("AAS53315")
            == "MFSTRICSLLARPFMVPIVPRFGSALLQKPLNGVVVPQFTRGFKVRTSVKKFCAHCYIVR"
            "RKGRVYVYCKSNNKHKQRQG"
        )
        assert (
            g.genetic_code.translate_dna(g.get_cds_sequence("AAS53315"))
            == "MFSTRICSLLARPFMVPIVPRFGSALLQKPLNGVVVPQFTRGFKVRTSVKKFCAHCYIVR"
            "RKGRVYVYCKSNNKHKQRQG"
        )

        assert (
            g.genetic_code.translate_dna(
                g.get_cds_sequence("AAS53315", g.df_proteins.loc["AAS53315"])
            )
            == "MFSTRICSLLARPFMVPIVPRFGSALLQKPLNGVVVPQFTRGFKVRTSVKKFCAHCYIVR"
            "RKGRVYVYCKSNNKHKQRQG"
        )
        with pytest.raises(ValueError):
            g.get_cds_sequence("AAS53315", g.df_proteins.loc["AAS53316"])

        assert (
            g.df_genes_meta.loc["AGOS_ADL186C"]["description"]
            == "Restriction of telomere capping protein 1 [Source:UniProtKB/Swiss-Prot;Acc:Q75AV6]"
        )

    def test_download_jobs_called_init(
        self, new_pipegraph, mock_download, shared_prebuild
    ):
        species = "Ashbya_gossypii"  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome(species, "41", prebuild_manager=shared_prebuild)
        g.find_prebuild("genome.fasta")  # this is the actual test.

    def test_species_formating(self, shared_prebuild):
        species = "ashbya_gossypii"  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        with pytest.raises(ValueError):
            EnsemblGenome(species, "41", prebuild_manager=shared_prebuild)

    def test_unknown_species_raises(self, mock_download, shared_prebuild):
        species = "Unknown_unknown"  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        EnsemblGenome(species, "41", prebuild_manager=shared_prebuild).download_genome()
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()

    def test_all_genes(self, mock_download, shared_prebuild):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()

        df = g.df_genes
        df2 = g.df_genes  # caching
        assert len(df) == 6910 + 4  # from the a2 locus
        assert df is df2
        assert df.loc["UMAG_12015"].strand == 1
        assert df.loc["UMAG_12015"].tss == 4370
        assert df.loc["UMAG_12015"].tes == 5366
        assert df.loc["UMAG_12015"].transcript_stable_ids == ("KIS66832",)

        assert df.loc["UMAG_00663"].strand == -1
        assert df.loc["UMAG_00663"].tss == 1_947_590
        assert df.loc["UMAG_00663"].tes == 1_945_040
        assert df.loc["UMAG_00663"].transcript_stable_ids == ("KIS72250",)

        assert df.loc["UMAG_03168"].transcript_stable_ids == ("KIS68597", "KIS68596")
        assert df["chr"].dtype.name == "category"
        assert df["biotype"].dtype.name == "category"

    def test_all_transcripts(self, mock_download, shared_prebuild):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()

        df = g.df_transcripts
        assert "gene_stable_id" in df.columns
        assert len(df) == 6928 + 4  # from the a2 locus
        assert df["chr"].dtype.name == "category"
        assert df["biotype"].dtype.name == "category"
        assert df.loc["KIS71021"].chr == "2"
        assert df.loc["KIS71021"].strand == 1
        assert df.loc["KIS71021"].start == 354_742
        assert df.loc["KIS71021"].stop == 356_690
        assert df.loc["KIS71021"].gene_stable_id == "UMAG_12118"
        assert df.loc["KIS71021"].biotype == "protein_coding"
        assert df.loc["KIS71021"].exons == ((354_742, 354_936), (355_222, 356_690))
        assert df.loc["KIS71021"].exon_stable_ids == ("KIS71021-1", "KIS71021-2")

    def test_df_exons(self, mock_download, shared_prebuild):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        exon_count = sum([len(x) for x in g.df_transcripts["exons"]])
        df_exons = g.df_exons
        assert len(df_exons) > 0
        assert len(g.df_exons) == exon_count
        assert hasattr(type(g).df_exons, "__call__")

    def test_all_proteins(self, mock_download, shared_prebuild):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()

        df = g.df_proteins
        assert df.strand.isin([1, -1]).all()

    def test_transcript_ids(self, mock_download, shared_prebuild):
        # test that ustilago has at least one gene with multilpe transcripts
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()

        df = g.df_genes
        assert (df.transcript_stable_ids.apply(lambda x: len(x)) > 1).any()

    def test_multiple_exon_transcripts(self, mock_download, shared_prebuild):
        # test that ustilago has at least one transcript with multiple exons
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()

        df = g.df_transcripts
        print(df.exons)
        assert (df.exons.apply(lambda x: len(x)) > 1).any()

    def test_get_additional_gene_gtfs(
        self, new_pipegraph, mock_download, shared_prebuild
    ):
        g = EnsemblGenome("Ashbya_gossypii", 33, shared_prebuild)
        assert len(g.get_additional_gene_gtfs()) == 0

        g = EnsemblGenome("Homo_sapiens", 74, shared_prebuild)
        assert "ribosomal_genes_grch37" in g.get_additional_gene_gtfs()[0].name
        assert g.get_additional_gene_gtfs()[0].exists()
        new_pipegraph.new_pipegraph()
        g = EnsemblGenome("Homo_sapiens", 75, shared_prebuild)
        assert "ribosomal_genes_grch38" in g.get_additional_gene_gtfs()[0].name
        assert g.get_additional_gene_gtfs()[0].exists()
        g = EnsemblGenome("Mus_musculus", 68, shared_prebuild)
        assert "ribosomal_genes_mm10" in g.get_additional_gene_gtfs()[0].name
        assert g.get_additional_gene_gtfs()[0].exists()
        new_pipegraph.new_pipegraph()
        g = EnsemblGenome("Mus_musculus", 67, shared_prebuild)
        assert len(g.get_additional_gene_gtfs()) == 0

    def test_get_additional_gene_gtfs_land_in_df_genes(
        self, mock_download, shared_prebuild
    ):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        assert len(g.get_additional_gene_gtfs()) == 1
        ppg.run_pipegraph()

        print(g.df_genes.chr.unique())
        assert "A2_pra2" in g.df_genes.index
        assert "A2_pra2.1" in g.df_transcripts.index

    def test_genes_iterator(self, mock_download, shared_prebuild):
        g = EnsemblGenome("Ashbya_gossypii", 41, shared_prebuild)
        ppg.run_pipegraph()
        genes = list(g.genes.values())
        assert len(genes) == len(g.df_genes)
        assert set([x.gene_stable_id for x in genes]) == set(g.df_genes.index)

    def test_transcript_iterator(self, mock_download, shared_prebuild):
        g = EnsemblGenome("Ashbya_gossypii", 41, shared_prebuild)
        ppg.run_pipegraph()
        transcripts = list(g.transcripts.values())
        assert len(transcripts) == len(g.df_transcripts)
        assert set([x.transcript_stable_id for x in transcripts]) == set(
            g.df_transcripts.index
        )

    def test_outside_of_ppg_after_download(self, mock_download, shared_prebuild):
        species = "Ashbya_gossypii"  # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome(species, "41", prebuild_manager=shared_prebuild)
        ppg.run_pipegraph()
        len_genes = len(g.df_genes)
        len_transcripts = len(g.df_transcripts)
        len_proteins = len(g.df_proteins)
        assert len_genes > 0
        assert len_transcripts > 0
        assert len_proteins > 0
        ppg.util.global_pipegraph = None
        g = EnsemblGenome(species, "41", prebuild_manager=shared_prebuild)
        assert len_genes == len(g.df_genes)
        assert len_transcripts == len(g.df_transcripts)
        assert len_proteins == len(g.df_proteins)

    def test_get_true_chromosomes(self, mock_download, shared_prebuild):
        # we need something with contigs and chromosomes
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        should = [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "23",
            "U37796.1",  # it is a true chromosome if there's chromosome: in the fasta description
        ]
        all_contigs = should + [
            "um_scaf_contig_1.256",
            "um_scaf_contig_1.265",
            "um_scaf_contig_1.271",
            "um_scaf_contig_1.264",
        ]
        assert set(g.get_true_chromosomes()) == set(should)
        assert set(g.get_chromosome_lengths()) == set(should + all_contigs)

    def test_get_true_chromosomes_genome_without_chromosomes(
        self, mock_download, shared_prebuild
    ):
        # we need something with contigs and chromosomes
        g = EnsemblGenome("Giardia_lamblia", 43, shared_prebuild)
        ppg.run_pipegraph()
        should = [
            "ctg02_1",
            "ctg02_2",
            "ctg02_3",
            "ctg02_4",
            "ctg02_5",
            "ctg02_6",
            "ctg02_7",
            "ctg02_8",
            "ctg02_9",
            "ctg02_10",
            "ctg02_11",
            "ctg02_12",
            "ctg02_13",
            "ctg02_14",
            "ctg02_15",
            "ctg02_16",
            "ctg02_17",
            "ctg02_18",
            "ctg02_19",
            "ctg02_20",
            "ctg02_21",
            "ctg02_22",
            "ctg02_23",
            "ctg02_24",
            "ctg02_25",
            "ctg02_26",
            "ctg02_27",
            "ctg02_28",
            "ctg02_29",
            "ctg02_30",
            "ctg02_31",
            "ctg02_32",
            "ctg02_33",
            "ctg02_34",
            "ctg02_35",
            "ctg02_36",
            "ctg02_37",
            "ctg02_38",
            "ctg02_39",
            "ctg02_40",
            "ctg02_41",
            "ctg02_42",
            "ctg02_43",
            "ctg02_44",
            "ctg02_45",
            "ctg02_46",
            "ctg02_47",
            "ctg02_48",
            "ctg02_49",
            "ctg02_50",
            "ctg02_51",
            "ctg02_52",
            "ctg02_53",
            "ctg02_54",
            "ctg02_55",
            "ctg02_56",
            "ctg02_57",
            "ctg02_58",
            "ctg02_59",
            "ctg02_60",
            "ctg02_61",
            "ctg02_62",
            "ctg02_63",
            "ctg02_64",
            "ctg02_65",
            "ctg02_66",
            "ctg02_67",
            "ctg02_68",
            "ctg02_69",
            "ctg02_70",
            "ctg02_71",
            "ctg02_72",
            "ctg02_73",
            "ctg02_74",
            "ctg02_75",
            "ctg02_76",
            "ctg02_77",
            "ctg02_78",
            "ctg02_79",
            "ctg02_80",
            "ctg02_81",
            "ctg02_82",
            "ctg02_83",
            "ctg02_84",
            "ctg02_85",
            "ctg02_86",
            "ctg02_87",
            "ctg02_88",
            "ctg02_89",
            "ctg02_90",
            "ctg02_91",
            "ctg02_92",
            "ctg02_93",
            "ctg02_94",
            "ctg02_95",
            "ctg02_96",
            "ctg02_97",
            "ctg02_99",
            "ctg02_98",
            "ctg02_100",
            "ctg02_101",
            "ctg02_102",
            "ctg02_103",
            "ctg02_104",
            "ctg02_105",
            "ctg02_106",
            "ctg02_107",
            "ctg02_108",
            "ctg02_109",
            "ctg02_110",
            "ctg02_111",
            "ctg02_112",
            "ctg02_113",
            "ctg02_114",
            "ctg02_115",
            "ctg02_116",
            "ctg02_117",
            "ctg02_118",
            "ctg02_120",
            "ctg02_119",
            "ctg02_121",
            "ctg02_122",
            "ctg02_123",
            "ctg02_124",
            "ctg02_125",
            "ctg02_127",
            "ctg02_126",
            "ctg02_129",
            "ctg02_128",
            "ctg02_130",
            "ctg02_131",
            "ctg02_132",
            "ctg02_133",
            "ctg02_134",
            "ctg02_135",
            "ctg02_136",
            "ctg02_137",
            "ctg02_138",
            "ctg02_139",
            "ctg02_140",
            "ctg02_141",
            "ctg02_142",
            "ctg02_144",
            "ctg02_143",
            "ctg02_145",
            "ctg02_146",
            "ctg02_147",
            "ctg02_148",
            "ctg02_149",
            "ctg02_151",
            "ctg02_150",
            "ctg02_152",
            "ctg02_153",
            "ctg02_155",
            "ctg02_154",
            "ctg02_157",
            "ctg02_156",
            "ctg02_158",
            "ctg02_159",
            "ctg02_160",
            "ctg02_161",
            "ctg02_162",
            "ctg02_163",
            "ctg02_164",
            "ctg02_165",
            "ctg02_167",
            "ctg02_168",
            "ctg02_166",
            "ctg02_169",
            "ctg02_170",
            "ctg02_171",
            "ctg02_172",
            "ctg02_174",
            "ctg02_176",
            "ctg02_175",
            "ctg02_178",
            "ctg02_177",
            "ctg02_179",
            "ctg02_180",
            "ctg02_181",
            "ctg02_183",
            "ctg02_182",
            "ctg02_184",
            "ctg02_185",
            "ctg02_186",
            "ctg02_187",
            "ctg02_188",
            "ctg02_189",
            "ctg02_190",
            "ctg02_191",
            "ctg02_193",
            "ctg02_192",
            "ctg02_194",
            "ctg02_196",
            "ctg02_195",
            "ctg02_197",
            "ctg02_198",
            "ctg02_200",
            "ctg02_199",
            "ctg02_201",
            "ctg02_202",
            "ctg02_203",
            "ctg02_204",
            "ctg02_206",
            "ctg02_208",
            "ctg02_207",
            "ctg02_205",
            "ctg02_209",
            "ctg02_210",
            "ctg02_211",
            "ctg02_212",
            "ctg02_214",
            "ctg02_213",
            "ctg02_215",
            "ctg02_216",
            "ctg02_220",
            "ctg02_218",
            "ctg02_221",
            "ctg02_217",
            "ctg02_219",
            "ctg02_223",
            "ctg02_222",
            "ctg02_224",
            "ctg02_225",
            "ctg02_226",
            "ctg02_227",
            "ctg02_228",
            "ctg02_229",
            "ctg02_231",
            "ctg02_230",
            "ctg02_232",
            "ctg02_234",
            "ctg02_235",
            "ctg02_233",
            "ctg02_237",
            "ctg02_238",
            "ctg02_236",
            "ctg02_239",
            "ctg02_241",
            "ctg02_240",
            "ctg02_242",
            "ctg02_245",
            "ctg02_243",
            "ctg02_246",
            "ctg02_244",
            "ctg02_247",
            "ctg02_249",
            "ctg02_250",
            "ctg02_248",
            "ctg02_252",
            "ctg02_251",
            "ctg02_253",
            "ctg02_254",
            "ctg02_255",
            "ctg02_256",
            "ctg02_257",
            "ctg02_258",
            "ctg02_259",
            "ctg02_260",
            "ctg02_173",
            "ctg02_261",
            "ctg02_262",
            "ctg02_263",
            "ctg02_266",
            "ctg02_265",
            "ctg02_268",
            "ctg02_264",
            "ctg02_267",
            "ctg02_269",
            "ctg02_271",
            "ctg02_270",
            "ctg02_273",
            "ctg02_272",
            "ctg02_274",
            "ctg02_275",
            "ctg02_276",
            "ctg02_277",
            "ctg02_278",
            "ctg02_279",
            "ctg02_280",
            "ctg02_282",
            "ctg02_281",
            "ctg02_283",
            "ctg02_284",
            "ctg02_285",
            "ctg02_286",
            "ctg02_287",
            "ctg02_289",
            "ctg02_288",
            "ctg02_290",
            "ctg02_291",
            "ctg02_292",
            "ctg02_293",
            "ctg02_294",
            "ctg02_295",
            "ctg02_296",
            "ctg02_297",
            "ctg02_298",
            "ctg02_299",
            "ctg02_300",
            "ctg02_301",
            "ctg02_302",
            "ctg02_303",
            "ctg02_304",
            "ctg02_305",
            "ctg02_306",
        ]

        assert set(g.get_true_chromosomes()) == set(should)
        assert set(g.get_true_chromosomes()) == set(g.get_chromosome_lengths())

    def test_newest_gene_ids(self, new_pipegraph, mock_download, shared_prebuild):
        # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        assert g.newest_stable_ids_for("UM05644") == set(["UMAG_05644"])
        assert g.newest_stable_ids_for("UMAG_05629") == set(["UMAG_05629"])
        assert g.newest_stable_ids_for("UM06501P0") == set([])
        assert g.newest_stable_ids_for("UM04933T0") == set([])
        with pytest.raises(KeyError):
            g.newest_stable_ids_for("no_such_gene")

    def test_get_external_dbs(self, new_pipegraph, mock_download, shared_prebuild):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        assert g.get_external_dbs() == [
            "BRENDA",
            "BioCyc",
            "ChEMBL",
            "ENA_FEATURE_GENE",
            "ENA_FEATURE_PROTEIN",
            "ENA_FEATURE_TRANSCRIPT",
            "ENA_GENE",
            "Ensembl_Fungi",
            "GO",
            "GOA",
            "IntAct",
            "IntEnz",
            "Interpro",
            "KEGG",
            "KEGG_Enzyme",
            "MEROPS",
            "MINT",
            "MetaCyc",
            "NCBI_TAXONOMY",
            "PHI",
            "PHIE",
            "PHIP",
            "PRIDE",
            "PUBMED",
            "PeroxiBase",
            "Reactome",
            "SWISS_MODEL",
            "UniParc",
            "UniPathway",
            "Uniprot/SPTREMBL",
            "Uniprot/SWISSPROT",
            "protein_id",
        ]

    def test_external_db_mapping(self, new_pipegraph, mock_download, shared_prebuild):
        # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        goa = g.get_external_db_to_gene_id_mapping("GOA")
        assert goa["A0A0D1CJ64"] == set(["UMAG_05734"])
        with pytest.raises(KeyError):
            g.get_external_db_to_gene_id_mapping("GOAnosuchthing")

    def test_external_db_mapping_transcript(
        self, new_pipegraph, mock_download, shared_prebuild
    ):
        # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        ena = g.get_external_db_to_gene_id_mapping("ENA_FEATURE_TRANSCRIPT")
        assert ena["CM003155.1:CDS:9690..12623"] == set(["UMAG_05624"])

    def test_external_db_mapping_translation(
        self, new_pipegraph, mock_download, shared_prebuild
    ):
        # the smallest eukaryotic species at the time of writing this at 2.8 mb
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        up = g.get_external_db_to_gene_id_mapping("Uniprot/SWISSPROT")
        assert up["P30598"] == set(["UMAG_10718"])
        upp = g.get_external_db_to_translation_id_mapping("Uniprot/SWISSPROT")
        assert upp["P30598"] == set(["KIS66849.N"])

    @patch("pypipegraph.util.checksum_file", return_value="5")
    def test_get_canonical_ids(self, new_pipegraph, mock_download, shared_prebuild):
        if hasattr(ppg, "is_ppg2"):
            import pypipegraph2 as ppg2

            g = EnsemblGenome("Homo_sapiens", 96, shared_prebuild, do_download=False)
            # g.job_genes()
            g._pb_find_server()
            g._pb_download_sql_table("gene")
            g._pb_download_sql_table("alt_allele")
            g._pb_download_sql_table("seq_region")
            g._pb_download_sql_table("seq_region_attrib")
            g._pb_download_sql_table("attrib_type")
            g._pb_download_gtf()
            g._pb_download_sql_table_definitions()
            g.job_genes()
            g.job_transcripts()
            ppg2.run()
            ppg2.change_global_pipegraph(None)
        else:
            g = EnsemblGenome("Homo_sapiens", 96, shared_prebuild)
            g._pb_find_server().callback()
            g._pb_download_sql_table("gene").callback()
            g._pb_download_sql_table("alt_allele").callback()
            g._pb_download_sql_table("seq_region").callback()
            g._pb_download_sql_table("seq_region_attrib").callback()
            g._pb_download_sql_table("attrib_type").callback()
            g._pb_download_gtf().callback()
            g._pb_download_sql_table_definitions().callback()
            g.job_genes().callback()
            g.job_transcripts().callback()
        assert g.name_to_canonical_id("DSEL") == "ENSG00000171451"
        assert g.name_to_canonical_id("THEMIS") == "ENSG00000172673"
        # test the breakage
        with pytest.raises(ValueError):
            g.name_to_canonical_id("SOD2")
        with pytest.raises(ValueError):
            g.name_to_canonical_id("IGF2")
        with pytest.raises(ValueError):
            g.name_to_canonical_id("ABCF2")
        with pytest.raises(ValueError):
            g.name_to_canonical_id("TBCE")

        assert g.name_to_canonical_id("SOD2", True) == "ENSG00000112096"
        assert g.name_to_canonical_id("IGF2", True) == "ENSG00000167244"
        assert g.name_to_canonical_id("ABCF2", True) == "ENSG00000033050"
        assert g.name_to_canonical_id("TBCE", True) == "ENSG00000284770"
        # assert g.name_to_canonical_id('HLA-DRB3') == 'ENSG00000230463'

    def test_same_same(self, new_pipegraph):
        g = EnsemblGenome("Ustilago_maydis", 33)
        g2 = EnsemblGenome("Ustilago_maydis", 33)
        assert g is g2
        new_pipegraph.new_pipegraph()
        g3 = EnsemblGenome("Ustilago_maydis", 33)
        assert g is not g3

    def test_additional_fasta(self, mock_download, shared_prebuild):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()
        seq = g.get_genome_sequence("U37796.1", 0, 100)
        assert (
            seq
            == (
                "taatcgtgaattgagctaggggcgccaagttacgtggcaaaagcgggctgactggcggcgaagatgtgt"
                "tggtctgcacctgagttcacgaacctgagac"
            ).upper()
        )

    def test_transcript_sequence(self, shared_prebuild):
        g = EnsemblGenome("Ustilago_maydis", 33, shared_prebuild)
        ppg.run_pipegraph()

        assert g.transcripts["KIS71708"].mrna == (
            "ATGAACGTCAAGCTTGCGCCGCCCGACGAGATGAACGGCGAGATCATTGCCATCCTCATT"
            "CTCTTCAACTTTCATTGCCGCGTTCAACTTAATGACGAGGCACAAGCCGCTCGCAAAAAG"
            "TTGCTCCTCTTCCTCATGGACAAAATTTACCAGACACGCGCACCTGCGCCTTCGTATGCA"
            "GCATTTGCAGACGAGCTCGGGGCCGCGATGGAAGGAGACGAGGACAATCACCGTATCATG"
            "ACCGACTACCTCGAGACCATGCTCGATCTGCTCCATGTTCCCGATGGTTTGACCAAGCTG"
            "TTCAACGAAAAGCTGAACAGAATCCTGCCTAGCTACGAGCCCATGGGATTGCTCAATGCC"
            "ACCGACATCTTCTTCGAGCGAAGGTCCTTCTTCGGCCTCTTTTTCCGTCGCATCAAGTTG"
            "ATCTTTGACAGTCTGGATTTGCAAACAAGGGATCACTTGACAATAGCTGCGCGCGCATGG"
            "AAAGAGGGACAGGCTTTCGACTTGAATGACAGCGAGCTCTCCGGTATCGACTCGGCGCAT"
            "TTGCTCGATGCAAGGCTTGGCGCGTTCCGAGACTACCAGCTTGGCTTACTCCGCGGCGAC"
            "TATACCATGGCAAAGGACAACATGGAAAGGTTCTTTGACTTTTACGCTCCTGGCGCTGAC"
            "CGCGAGCTTCATCAACACACGCTTCTTCATCTAGCCGCTTTCCATGTCAGGACTGAGAGC"
            "TTCTGCGCAGCCAAAGCAGCCTTGGACGAGGCTATAAGCCTTGCGCGGTCGGCAAACGAC"
            "AGCGAGTGCATCTCTGCCTGCGAGAGCCTCATGCAATACATCCAAGGTGTGGGCACTAGC"
            "ACCCTCGCTTCTGTTCCGGGAGCAGCCAATGTCTCGACGAACGAGAGGCAACGCCGACCT"
            "ATCTATGACGCGGTGTGGCAGACCCGTTGCAGCCTAGCTAAGGGCCGCTCAGCAATCGAG"
            "ACACTTCAAGACTTGGAGGATTGTGCTGCACCTTCTCAGCCTTCTAGGGACTCGTTGGCG"
            "GCCAGCGAGGCCTCCTTGCAGCTCATTCAGGATGCCAAACGACGTCTCGGCAGGGATACA"
            "TTTCCAAGCGACGGTGAAGTCGCGCGTCTTTGGGATACCTTGGGCCAGCCAGCGCTTGCT"
            "GACGTCTATCGGAATCGGAACATTGCAGGTGCAAATGGTCGAGCACGCTCAGCTTTGCAA"
            "GAAGAAAGCCGAATCGATTGCATCTGTCACAAAGCCAAAACG"
            "CTTGCACGTGCGGGAGAGTACGAAGCGGCGCTCAGCTTGCTGGTCTCACCGGCCACTTTC"
            "GAAGCGATCTCATTCTGCGAGTACACGATCTGGCACAGGGCAATCGCTGAGGTCCTGCGC"
            "TTACGAGCAACGCGCAGACAAGATGTAGCGACCTTGCAATTGTTGGCAGAGAGCCTGCCG"
            "GGTTCCGACCAGGCTCATGTCGATCGCGACGTTGAAGATGCTGTCGATTCGCCGAGCGCT"
            "TTGGTCGAGCTCGCGTTGCGCTGCCTCAAGTCCGGAAAATCGAGTGCGACTGAGAAGAAT"
            "TTCGGATTGCGATTCAAAGAGTACAGCCGACTTACAGCAGAGCAAGTGGCAGAAGCGCTG"
            "CTGAACAAAGCTGCGATAAGGATGAAACGTGGTCGACCGATCCTCGCACTGATGCCTACA"
            "CTGGCAAGTCTATCGATTGCCAAGGACATGGAGTGCAATCGGTTGATCCTCAATGCAAGA"
            "GTTCAGCTAGCTGAGGCATTGGGGCTTCAGCTCAAGATGCCAGATGGAGCACGGCTGCTC"
            "TTGGAATCGGATCTGCCCAACTGCCTGTCGAGTGATGATGTCGAGCTGAGAGCGCGCGCG"
            "AAGTGGACGTATGCACGGATGCTGCTCTCATGCTCGGACAAGCAAGAGCGCGAAGATCTG"
            "ACCAAGGTGCTCTACTGGTTGCGAGAAGCCGAAAGAG"
            "ACGCACAACAGGCCGAATGCCTCGAGCTTCACACGCAGATTCTCTACTATATGTTGCGGC"
            "TGCACCACCACCTGGGCGATGACAGAGAAACAATCTCTGTGACAGCTCGTTTAGATACAG"
            "TAGAGCGCGCTTGGACTCGCTTGGATGCTTCGCAGGATCAAGCCCATCTGCAACAGGTTC"
            "GCCAAATCCTAGATATCGTTGTCTCTGTCGCAGGTTATGTGGCTAGCGGAGAAGCTGCGA"
            "ACAAGCGTCTTGAAATGGTCTAG"
        )

        assert (
            g.transcripts["KIS71709"].mrna == "ATGGCATTTTCAGAAGATACCAAG"
            "GAGCGAATCATTAAGGCGGTCGATGTCTCCAAGACTTTGTTGCATTACGGCTG"
            "GGTGCCTTTCGTTCTTTACATCGGCTTCACCCGAAGCACGCCCCAGCCTAGCTTGATCAA"
            "GCTCATCAGTCCTCTCGCATGA"
        )
        assert g.transcripts["KIS71709"].mrna == g.transcripts["KIS71709"].cdna
