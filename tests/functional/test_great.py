import mbf.genomes
import collections
import pickle
import os
import pandas as pd
import pytest
import mbf.functional.great
import gzip
import pypipegraph2 as ppg
from pathlib import Path
from mbf.genomes import HardCodedGenome

data_path = Path(__file__).parent


default_chr_lengths = {
    "1": 100_000,
    "2": 200_000,
    "3": 300_000,
    "4": 400_000,
    "5": 500_000,
}


def DummyGenome(df_genes, df_transcripts=None, chr_lengths=default_chr_lengths):
    assembly = "hg19"

    df_genes = df_genes.rename(columns={"stable_id": "gene_stable_id"})
    if not "start" in df_genes.columns:
        starts = []
        stops = []
        for idx, row in df_genes.iterrows():
            if row["strand"] == 1:
                starts.append(row["tss"])
                stops.append(row["tes"])
            else:
                starts.append(row["tes"])
                stops.append(row["tss"])
        df_genes = df_genes.assign(start=starts, stop=stops)
    if not "biotype" in df_genes.columns:
        df_genes = df_genes.assign(biotype="protein_coding")
    if not "name" in df_genes.columns:
        df_genes = df_genes.assign(name=df_genes.index)
    df_genes = df_genes.sort_values(["chr", "start"])
    df_genes = df_genes.set_index("gene_stable_id")
    if df_transcripts is None:
        df_transcripts = df_genes.copy().reset_index()
        df_transcripts = df_transcripts.assign(
            transcript_stable_id=["TR" + x for x in df_transcripts["gene_stable_id"]]
        )

    if not "biotype" in df_transcripts.columns:
        df_transcripts = df_transcripts.assign(biotype="protein_coding")
    if not "name" in df_transcripts.columns:
        df_transcripts = df_transcripts.assign(name=df_transcripts.index)

    if "exons" not in df_transcripts.columns:
        df_transcripts = df_transcripts.assign(
            exons=[
                [(start, stop)]
                for start, stop in zip(df_transcripts["start"], df_transcripts["stop"])
            ]
        )

    if "exons" in df_transcripts.columns:
        if len(df_transcripts["exons"].iloc[0]) == 3:
            df_transcripts = df_transcripts.assign(
                exons=[(x[0], x[1]) for x in df_transcripts["exons"]]
            )
        exon_stable_ids = []
        for row in df_transcripts["exons"]:
            here = [tuple(["exon_%s_%i" % (x, ii)]) for (ii, x) in enumerate(row)]
            exon_stable_ids.append(here)
        df_transcripts = df_transcripts.assign(exon_stable_ids=exon_stable_ids)
    df_transcripts = df_transcripts.set_index("transcript_stable_id")
    res = HardCodedGenome("dummy", chr_lengths, df_genes, df_transcripts, None)
    res.assembly = assembly
    return res


class ComparisonFunctionalGroups:

    """The comparison groups from Cory McLean"""

    name = "comparison_groups"

    def get_sets(self, genome):
        op = gzip.GzipFile(data_path / "great_data" / "hg18.regdoms.gz", "rb")
        s = set()
        for line in op:
            line = line.decode("utf-8")
            line = line.strip()
            if line:
                line = line.split()
                stable_id = line[-1]
                s.add(stable_id)
        return {'GO:0015629 ("actin cytoskeleton")': s}

    def get_url(self, set):
        return ""


class FakeGenomeForComparison:

    """A genome that has just the genes Cory McLean provided for us"""

    name = "compgenome"

    def download_genome(self):
        return []

    def get_chromosome_lengths(self):
        bed = mbf.fileformats.bed.read_bed(
            data_path / "great_data" / "hg18.antigap.bed.gz"
        )
        max_by_chr = collections.defaultdict(int)
        for entry in bed:
            # max_by_chr[entry.refseq] = max(max_by_chr[entry.refseq], entry.position + entry.length)
            max_by_chr[entry.refseq] += entry.length
        return dict(max_by_chr.items())  # not a defaultdict...

    def get_all_genes(self, filter_non_canonical_chromosomes=False):
        op = gzip.GzipFile(os.path.join(data_path, "great_data", "hg18.loci.gz"), "r")
        data = {
            "chr": [],
            "tss": [],
            "tes": [],
            "strand": [],
            "stable_id": [],
            "name": [],
        }
        for line in op:
            line = line.strip().decode("utf-8")
            if line:
                line = line.split()
                data["stable_id"].append(line[0])
                data["name"].append(line[-1])
                data["chr"].append(line[1])
                data["tss"].append(int(line[2]))
                if line[3] == "-":
                    data["strand"].append(-1)
                    data["tes"].append(int(line[2]) - 1)
                elif line[3] == "+":
                    data["strand"].append(1)
                    data["tes"].append(int(line[2]) + 1)
                else:
                    raise ValueError("not a valid strand: %s" % line)

        df = pd.DataFrame(data)
        assert df["stable_id"].dtype == "object"
        assert df["name"].dtype == "object"
        assert isinstance(df["name"].iloc[0], str)
        assert isinstance(df["stable_id"].iloc[0], str)
        df.index = df["stable_id"]
        return df

    @property
    def genes(self):
        res = self.get_all_genes(True)
        res = {x["stable_id"]: x for _, x in res.iterrows()}
        return res

    def get_dependencies(self):
        return []

    def gene_id_to_name(self, stable_id):
        return [stable_id]


def Fake_GREAT_RegulatorRegion_(fake_genome):
    """Fake regulatory region from Cory McLean bed file..."""

    def load_data():
        df = pd.read_csv(
            data_path / "great_data" / "hg18.regdoms.gz", sep="\t", header=None
        )
        df.columns = [
            "chr",
            "start",
            "stop",
            "name",
        ]
        df = df.assign(name=df["name"].astype(str))
        return df

    return mbf.genomics.regions.GenomicRegions(
        "GREAT_RegulatorRegion_Mc_Lean", load_data, [], fake_genome, on_overlap="ignore"
    )


class TestGreat:
    def test_broken_regulatory_region_building(self, new_pipegraph):
        genome = DummyGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            )
        )
        # g = mbf.genomics.genes.Genes(genome)
        # now, even though we do gap handling, this should not change the
        # generated regions
        regions = mbf.functional.great.GreatRegulatoryRegions(genome)
        regions.load()

        def dump(of):
            # since we need this to be loaded to test.
            print(regions.regulatory_region_by_gene)

            # this is what I wrote after reading the paper.
            # but it's not what their createRegulatoryDomains c code actually does.
            assert regions.regulatory_region_by_gene["fake1"] == [("1", 0, 6000)]
            assert regions.regulatory_region_by_gene["fake2"] == [("1", 4400, 100000)]
            assert regions.regulatory_region_by_gene["fake3"] == [("2", 0, 200000)]
            of.write_text("sentine")

        ppg.FileGeneratingJob("force_load", dump).depends_on(
            mbf.functional.great._add_regulatory_regions_by_gene_to_gr(regions)
        )
        # regions.dump("cache/shu.tsv")
        # genome.get_non_gaps().load()
        ppg.run()

    @pytest.mark.skip("This test is not implemented")
    def test_gap_handling_works(self, new_pipegraph):
        raise NotImplementedError("This test is not implemented")
        genome = DummyGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                    },
                ]
            ),
            gaps=pd.DataFrame(
                [
                    {"chr": "1", "start": 0, "stop": 500},
                    {"chr": "2", "start": 0, "stop": 500},
                    {"chr": "2", "start": 4000, "stop": 4567},
                ]
            ),
        )
        # g = mbf.genes.Genes(genome)
        # now, even though we do gap handling, this should not change the
        # generated regions
        regions = mbf.functional.great.GREAT_RegulatorRegion(genome, exclude_gaps=True)
        regions.load()
        # regions.dump("cache/shu.tsv")
        genome.get_non_gaps().load()
        ppg.run()
        assert regions.regulatory_region_by_gene["fake1"] == [("1", 500, 6000)]
        assert regions.regulatory_region_by_gene["fake2"] == [("1", 4400, 100000)]
        assert regions.regulatory_region_by_gene["fake3"] == [
            ("2", 500, 4000),
            ("2", 4567, 200000),
        ]

    # this test takes about 30 second.s..
    def test_regulatory_region_calculation_5kb_1k_basal_1mb_extension(
        self, new_pipegraph
    ):
        """test_regulatory_region_calculation_5kb_1k_basal_1mb_extension e compare our calculation on a genome that was also run against GREAT's createRegulatoryDomains previously"""
        import time

        print("start", time.time())
        chr_lengths = {
            "1": 249250621,
            "10": 135534747,
            "11": 135006516,
            "12": 133851895,
            "13": 115169878,
            "14": 107349540,
            "15": 102531392,
            "16": 90354753,
            "17": 81195210,
            "18": 78077248,
            "19": 59128983,
            "2": 243199373,
            "20": 63025520,
            "21": 48129895,
            "22": 51304566,
            "3": 198022430,
            "4": 191154276,
            "5": 180915260,
            "6": 171115067,
            "7": 159138663,
            "8": 146364022,
            "9": 141213431,
            "MT": 16569,
            "X": 155270560,
            "Y": 59373566,
        }
        df = pd.read_csv(
            os.path.join(
                data_path, "great_data", "ensembl_homo_sapiens_60_genes.tsv.gz"
            ),
            sep="\t",
        )
        df = df.assign(
            start=[
                row["tss"] if row["tss"] < row["tes"] else row["tes"]
                for dummy_idx, row in df.iterrows()
            ],
            stop=[
                row["tes"] if row["tss"] < row["tes"] else row["tss"]
                for dummy_idx, row in df.iterrows()
            ],
        )
        df["chr"] = df["chr"].astype("str")
        # df = df[df.where(lambda row: row['chr'] in chr_lengths), :]
        df = df[df["chr"].isin(chr_lengths.keys())]
        genome = DummyGenome(df, chr_lengths=chr_lengths)
        # since the C version does not exclude them either...
        reg_regions = mbf.functional.great.GreatRegulatoryRegions(genome)
        # reg_regions.dump('results/regregions.bed')
        # reg_regions.regions.build_intervals()
        reg_regions.load()

        def dump(of):
            import pickle

            pickle.dump(reg_regions.regulatory_region_by_gene, open(of, "wb"))

        jj = mbf.functional.great._add_regulatory_regions_by_gene_to_gr(reg_regions)
        jx = ppg.FileGeneratingJob("force_load", dump).depends_on(jj)
        print("run pipeline", time.time())
        ppg.run()
        # hack around ppg deleting the attribute after run
        reg_regions.regulatory_region_by_gene = pickle.load(open(jx.job_id, "rb"))
        print("comp", time.time())
        op = gzip.GzipFile(
            data_path
            / "great_data"
            / "c_great_reg_regions_for_homo_sapiens_ensembl_60.bed.gz",
            "r",
        )
        # first we check whether we created the very same regulatory regions...
        # note: Non of these are split in two by the gaps...
        print(reg_regions.regulatory_region_by_gene["ENSG00000223972"])
        any_failed = False
        for line in op:
            line = line.decode("utf-8")
            if line:
                line = line.strip().split()
                chr, start, stop, gene, dummy, strand = line
                start = int(start)
                stop = int(stop)
                supposed = [(chr, start, stop)]
                actual = reg_regions.regulatory_region_by_gene[gene]
                if supposed != actual:
                    print(supposed)
                    print(actual)
                    print(' ')
                    any_failed = True
                # print(len(reg_regions.get_overlapping(chr, start, stop)))
                assert (
                    gene in reg_regions.get_overlapping(chr, start, stop)["name"].values
                )
        op.close()
        if any_failed:
            raise ValueError("discrepancies")
        print("done", time.time())

    def test_against_cory_mclean_data(self, new_pipegraph):
        """
        Test against data from Cory McLean, the original author of great who has gratiously provided us with some sample data
        so that we may compare to this implementation"""
        genome = FakeGenomeForComparison()
        query = mbf.genomics.regions.GenomicRegions_FromBed(
            "SRF", os.path.join(data_path, "great_data", "SRF.sigPeaks.bed.gz"), genome
        )
        fg = ComparisonFunctionalGroups()
        rregions = Fake_GREAT_RegulatorRegion_(genome)
        # mbf.functional.great._GREAT_RegulatorRegion_cache[
        #     genome.name + "basal_5k+1kb_extension_1MB" + "exclude"
        # ] = rregions
        # now, we say 'exclude gaps', but in truth that never happens in this
        # dataset...'
        antigap = mbf.genomics.regions.GenomicRegions_FromBed(
            "Fake_non_gaps", data_path / "great_data" / "hg18.antigap.bed.gz", genome
        )

        great = mbf.functional.great.GREAT(
            query_gis=query,
            function_gene_groups_or_list_of_such=fg,
            regulatory_regions=rregions,
            non_gap_regions=antigap,
        )
        great.write("results/great.tsv")
        ppg.run()

        my_result = pd.read_csv("results/great.tsv", sep="\t")
        their_result = pd.read_csv(
            data_path
            / "great_data"
            / "20101208_GREAT_web_output_all-GOCellularComponent.xls.gz",
            sep="\t",
        )
        their_result = their_result.set_index("GO TERM")
        # self.assertAlmostEqual(
        assert (
            my_result.iloc[0]["p-value binomial"]
            - their_result.loc["GO:0015629"]["Binom Raw P-Value"]
        ) < 1e-10
