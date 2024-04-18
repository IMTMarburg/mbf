import pypipegraph2 as ppg
import scipy.stats
import re
import collections
import hashlib
import pandas as pd
import subprocess
import mbf.externals
import pandas
from statsmodels.stats.multitest import multipletests
from pathlib import Path
from mbf.genomics.regions import GenomicRegions_FromBed


"""Great is a 'statistical' test for gene set enrichment of chipseq peaks.

It turns a genome into regulatory regions.
And it takes the anti gap regions.

Then for eache annotated gene set:
    Numerator = (UNION of all regulatory domains of genes annotated with the term) INTERSECT (UNION of all antigap ranges).
    Demonitator = is UNION of all anti gap ranges

"""


def _prepare_great_regulatory_regions(
    genome,
    params,
):
    """Actually call the great binary to create the regulatory regions.

    See GreatRegulatoryRegions for details on the parameters.
    """

    allowed_methods = ["basalPlusExtension", "oneClosest", "twoClosest"]
    if not params["method"] in allowed_methods:
        raise ValueError("method must be one of {}".format(allowed_methods))
    cache_file = Path(
        f"cache/GREAT/{genome.name}_{params['method']}_{params['maxExtension']}_{params['basalUpstream']}_{params['basalDownstream']}.bed"
    )

    try:
        nix_path = mbf.externals.externals.get_nix_store_path_from_binary(
            "createRegulatoryDomains"
        )

        great_binary = nix_path / "bin" / "createRegulatoryDomains"
        if hasattr(genome, "assembly"):
            great_reg_file = (
                nix_path / f"curated_regulatory_regions/{genome.assembly}.tsv"
            )
            if not great_reg_file.exists:
                great_reg_file = None
        else:
            great_reg_file = None

    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(
            "createRegulatoryDomains not found in path. Please install GREAT (IMTMarburg/flakes/GREAT)",
            e,
        )

    def gen(cache_file):
        try:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            tf_input = cache_file.with_name(cache_file.name + "_input.tss")
            tf_chrom_sizes = cache_file.with_name(cache_file.name + "chrom.sizes")
            tf_chrom_sizes.write_text(
                "\n".join(
                    f"{k}\t{v}" for k, v in genome.get_chromosome_lengths().items()
                )
            )
            # TSS.in
            # This is a file holding a list of all genes to which you want to assign regulatory domains.  Each line of the file
            # should correspond to a single gene to which you assign a regulatory domain.  Each line should have four fields
            # tab-delimited:
            #           chromosome      transcription start site      strand      geneName

            # so... what about multiple TSS?
            # -> Note that GREAT measures all distances from the transcription start site of a gene's canonical isoform.
            # and how to integerate the
            # what is strand, plus/minus 1, + or -?? Source says + or -
            if great_reg_file:
                df = pd.read_csv(
                    great_reg_file, sep="\t", header=None, index_col=False, comment="#"
                )
                df[0] = df[0].str[3:]  # strip chr
                lookup = {row[3]: row for row in df.itertuples(index=False)}
            else:
                lookup = {}

            try:
                gtf = genome.get_gtf(["transcript"])["transcript"]  # we need the tags.
                canonical_transcripts = set(
                    gtf["transcript_id"][
                        (gtf["tag0"] == "Ensembl_canonical")
                        | (gtf["tag1"] == "Ensembl_canonical")
                        | (gtf["tag2"] == "Ensembl_canonical")
                        | (gtf["tag3"] == "Ensembl_canonical")
                        | (gtf["tag4"] == "Ensembl_canonical")
                        | (gtf["tag5"] == "Ensembl_canonical")
                    ]
                )
            except TypeError:
                # non ensembl genomes -> everything
                canonical_transcripts = set(genome.transcripts.keys())

            # now for each transcript, if it's canonical...
            with open(tf_input, "w") as input_file:
                for transcript in genome.transcripts.values():
                    if transcript.transcript_stable_id in canonical_transcripts:
                        # get the gene name
                        gene_stable_id = transcript.gene_stable_id
                        chrom = transcript.chr
                        start = transcript.start
                        stop = transcript.stop
                        if transcript.strand not in (1, -1):
                            raise ValueError("Strand unexpected", transcript)
                        strand = "+" if transcript.strand == 1 else "-"
                        if strand == "+":
                            tss = start
                        else:
                            tss = stop
                        input_file.write(
                            f"{chrom}\t{tss}\t{strand}\t{gene_stable_id}\n"
                        )
                        # print(f"{chrom}\t{tss}\t{strand}\t{gene_stable_id}\n")
            tf = cache_file.with_suffix(".temp")
            subprocess.check_call(
                [
                    "createRegulatoryDomains",
                    tf_input.absolute(),
                    tf_chrom_sizes.absolute(),
                    params["method"],
                    f'-maxExtension={params["maxExtension"]}',
                    f'-basalUpstream={params["basalUpstream"]}',
                    f'-basalDownstream={params["basalDownstream"]}',
                    tf.absolute(),
                ]
            )
            tf2 = cache_file.with_suffix(".temp2")
            with open(tf2, "w") as output:
                with open(tf, "r") as input:
                    for line in input.readlines():
                        chrom, start, stop, gene, tss, strand = line.split("\t")
                        if gene in lookup:
                            strand = genome.genes[gene].strand
                            strand = "+" if strand == 1 else "-"
                            chr, start, stop, _ignored_gene, tss = lookup[gene]
                            output.write(
                                "\t".join(
                                    [
                                        str(x)
                                        for x in [chr, start, stop, gene, tss, strand]
                                    ]
                                )
                                + "\n"
                            )
                        else:
                            output.write(f"{line}")
            tf2.rename(cache_file)
        finally:
            if "tf" in locals():
                if tf.exists():
                    tf.unlink()
            if "tf_input" in locals():
                if tf_input.exists():
                    tf_input.unlink()
            if "tf_chr_sizes" in locals():
                if tf_chrom_sizes.exists():
                    tf_chrom_sizes.unlink()
            if "tf2" in locals():
                if tf2.exists():
                    tf2.unlink()

    job = (
        ppg.FileGeneratingJob(cache_file, gen)
        .depends_on(genome.download())
        .depends_on_params({"binary": great_binary})
        .self
    )
    if great_reg_file:
        job.depends_on_file(great_reg_file)
    return job


_great_regulatory_regions = {}


def GreatRegulatoryRegions(
    genome,
    params={
        "maxExtension": 1000000,
        "basalUpstream": 5000,
        "basalDownstream": 1000,
        "method": "basalPlusExtension",
    },
    exclude_gaps=False,
):
    name = f"great__{genome.name}_{params['method']}_{params['maxExtension']}_{params['basalUpstream']}_{params['basalDownstream']}"
    if not name in _great_regulatory_regions:
        prep_job = _prepare_great_regulatory_regions(genome=genome, params=params)
        # we'ell use this one as input, for it is sorted (which might not be true for the patched version we created above)
        res = GenomicRegions_FromBed(
            name, prep_job, genome, summit_annotator=False, on_overlap="ignore"
        )
        if exclude_gaps:
            raise ValueError("Todo")
        _great_regulatory_regions[name] = res
    return _great_regulatory_regions[name]


def check_function_gene_groups_or_list_of_such(function_gene_groups_or_list_of_such):
    """Some functional methods take either a GeneGroups(interface) object, or a list of such.
    This makes sure they all have the right functions to be used.
    """
    if hasattr(function_gene_groups_or_list_of_such, "get_sets"):
        result = [function_gene_groups_or_list_of_such]
        failed = False
    elif isinstance(function_gene_groups_or_list_of_such, list):
        result = function_gene_groups_or_list_of_such
        failed = False
        for f in function_gene_groups_or_list_of_such:
            if not hasattr(f, "get_sets"):
                failed = True
    else:
        failed = True
    if failed:
        raise ValueError(
            "function_gene_groups_or_list_of_such must be either something conforming to the FunctionalGeneGroups interface, or a list (not a collection) of such objects"
        )
    return result


_gap_regions = {}


def calc_gap_regions(genome):
    if genome not in _gap_regions:

        def find_gaps():
            result = collections.defaultdict(list)
            for name, length in genome.get_chromosome_lengths().items():
                seq = genome.get_genome_sequence(name, 0, length)
                for match in re.finditer("N{1000,}", seq):
                    result["chr"].append(name)
                    result["start"].append(match.start(0))
                    result["stop"].append(match.end(0))
            return pd.DataFrame(result)

        inner = mbf.genomics.regions.GenomicRegions(
            f"Gaps_{genome.name}_calc",
            find_gaps,
            [],
            genome,
        ).write_bed()[
            0
        ]  # job only
        outer = mbf.genomics.regions.GenomicRegions_FromBed(
            f"Gaps_{genome.name}",
            inner,
            genome,
            summit_annotator=False,
            on_overlap="raise",
        )
        _gap_regions[genome] = outer
    return _gap_regions[genome]


_non_gap_regions = {}


def calc_non_gap_regions(genome):
    if not genome in _non_gap_regions:
        gaps = calc_gap_regions(genome)
        r = gaps.convert(
            f"NonGaps_{genome.name}", mbf.genomics.regions.convert.invert()
        )
        r.write_bed()  # todo: is this necessary?
        _non_gap_regions[genome] = r
    return _non_gap_regions[genome]


def _add_regulatory_regions_by_gene_to_gr(regulatory_regions):
    def load():
        regulatory_region_by_gene = {}
        for dummy_idx, row in regulatory_regions.df.iterrows():
            if not row["name"] in regulatory_region_by_gene:
                regulatory_region_by_gene[row["name"]] = []
            regulatory_region_by_gene[row["name"]].append(
                (row["chr"], row["start"], row["stop"])
            )
        return regulatory_region_by_gene

    return ppg.AttributeLoadingJob(
        regulatory_regions.name + "_fill_regulatory_region_by_gene",
        regulatory_regions,
        "regulatory_region_by_gene",
        load,
    ).depends_on(regulatory_regions.load())


def hypergeom_probability_larger(x, r, b, n):
    """Return p(X >= x) on a hypergeometric distribution
    x = number of white balls drawn,
    r = number of white balls in the urn
    b = number of black balls in the urn
    n = number of balls drawn.
    Or as contingency diagram:
        white drawn (x) | black drawn     | # drawn(n)
        white not drawn | black not drawn | # not drawn
        -------------------------------------------
        total white (r) | total black (b) | grand total


    """
    return scipy.stats.hypergeom(
        M=r + b,  # total number of balls
        n=r,  # number of white balls
        N=n,  # no of balls drawn
    ).sf(
        x - 1  # no of white balls drawn
    )


def fdr_control_benjamini_hochberg(
    df, p_value_column, fdr_output_column, drop_output_if_exist=False
):
    reject, pvals_corrected, _, _ = multipletests(
        df[p_value_column], method="fdr_bh", alpha=0.05
    )
    if fdr_output_column in df.columns:
        if drop_output_if_exist:
            df = df.drop(columns=[fdr_output_column])
        else:
            raise ValueError(f"Column {fdr_output_column} already exists in dataframe")
    return df.assign(**{fdr_output_column: pvals_corrected})


class GREAT:
    """GREAT implementation

    Follows McLean et al, "Great improves functional interpretation of cis-regulatory regions,
    Nature Biotechnology, Volume 28, 5, May 2010 doi:10.1038/nbt.1630.

    Basically, it defines regulatory regions for each gene, and scores binding sites(not gene sets!)
    within the regions of one particular group of genes by either a binomial (considering the size
    of the regulatory regions) and hypergeometric (which ignores the prior imposed by the size of these regions) test
    """

    def __init__(
        self,
        query_gis,
        function_gene_groups_or_list_of_such,
        regulatory_regions=None,
        non_gap_regions=None,
    ):
        self.query_gis = query_gis
        self.list_of_gene_groups = check_function_gene_groups_or_list_of_such(
            function_gene_groups_or_list_of_such
        )
        self.name = "great_%s_vs_%s" % (
            self.query_gis.name,
            ",".join(x.name for x in self.list_of_gene_groups),
        )
        ppg.assert_uniqueness_of_object(self)
        self.cache_path = Path("cache", "functional", "GREAT")
        self.cache_path.mkdir(parents=True, exist_ok=True)
        self.non_gap_regions = (
            calc_non_gap_regions(self.query_gis.genome)
            if non_gap_regions is None
            else non_gap_regions
        )
        if regulatory_regions is None:
            self.regulatory_regions = GreatRegulatoryRegions(self.query_gis.genome)
        else:
            self.regulatory_regions = regulatory_regions

    def perform(self):
        """Run GREAT analysis (job)"""

        def calc():
            data = collections.defaultdict(list)
            # technially, we should remove the gaps from each one...
            # genome_base_count = sum(
            #     [
            #         int(x)
            #         for x in self.query_gis.genome.get_chromosome_lengths().values()
            #     ]
            # )
            genome_base_count = self.non_gap_regions.covered_bases

            query_hit_genes = self.find_genes_hit_by_query()
            all_genes = set(self.query_gis.genome.genes.keys())
            print("in testing", self.query_gis.name)
            print(
                "The test set of %i genomic regions picked %i genes"
                % (self.query_gis.get_no_of_entries(), len(query_hit_genes))
            )
            print("len all_genes", len(all_genes))
            print(query_hit_genes)
            assert isinstance(list(all_genes)[0], str)
            assert self.list_of_gene_groups
            for group in self.list_of_gene_groups:
                sets = group.get_sets(self.query_gis.genome)
                all_in_sets = set()
                for setname, genes in sets.items():
                    all_in_sets.update(genes)
                # print('%s has %i terms covering %i genes' % (group.name, len(sets), len(all_in_sets))
                for setname, genes in sets.items():
                    assert isinstance(list(genes)[0], str)
                    try:
                        # filter out those that the set has, but that are
                        # not(no longer) in the genome.
                        prev_genes = genes
                        genes = set(genes).intersection(all_genes)
                        if len(genes) == 0:
                            print("prev", prev_genes)
                            print("all", all_genes)
                            print("filtered all genes from set", setname)
                    except:  # noqa: E722
                        print(setname)
                        raise
                    set_hits = 0
                    set_hit_genes = set()
                    intervals_seen = set()
                    for gene, interval_id in query_hit_genes:
                        if (
                            interval_id in intervals_seen
                        ):  # GREAT apperently only counts each interval once within a geneset, which I guess makes sense, or it would introduce a bias to gene-dense-closely-annotated clusters
                            continue
                        intervals_seen.add(interval_id)
                        if gene in genes:
                            set_hits += 1
                            set_hit_genes.add(gene)
                    hypergeom_hit_genes = set([x[0] for x in query_hit_genes])
                    overlap = hypergeom_hit_genes.intersection(genes)
                    covered_bases = self.calculate_regulator_region_size(genes)
                    n = self.query_gis.get_no_of_entries()
                    p = float(covered_bases) / genome_base_count
                    k = set_hits
                    p_value_binomial = scipy.stats.binom.sf(k - 1, n, p)
                    # print(setname
                    # print('print(white pulled', len(overlap)
                    # print('white in urn', len(genes)
                    # print('black in urn',len(self.regulatory_regions.regulatory_size) - len(genes)
                    # print('drawn',n, 'hit', len(query_hit_genes)
                    p_value_hypergeom = hypergeom_probability_larger(
                        len(overlap),  # white balls drawn
                        len(genes),  # white balls in urn
                        len(all_genes) - len(genes),  # black balls in urn
                        len(hypergeom_hit_genes),  # balls drawn
                    )
                    # print('p_hypergom', p_value_hypergeom
                    # print(''
                    data["Group"].append(group.name)
                    data["Set"].append(setname)
                    data["Set Size"].append(len(genes))
                    data["Set hypergeom hits"].append(len(overlap))
                    data["hypergeom white balls"].append(len(genes))
                    data["hypergeom black balls"].append(len(all_genes) - len(genes))
                    data["hypergeom drawn balls"].append(len(hypergeom_hit_genes))
                    data["Hits"].append(k)
                    data["Drawn"].append(n)
                    data["p-value binomial"].append(p_value_binomial)
                    data["p-value hypergeometric"].append(p_value_hypergeom)
                    data["set covered bases"].append(covered_bases)
                    data["genome bases"].append(genome_base_count)
                    data["p"].append(p)
                    data["hit genes"].append(
                        ", ".join(
                            sorted(
                                self.query_gis.genome.genes[gene_stable_id].name
                                for gene_stable_id in set_hit_genes
                            )
                        )[:16000]
                    )
                    data["Link"].append(group.get_url(setname))
                    # data['p-value great direct'].append(self.call_great(genes, n, k))
            df = pd.DataFrame(data)
            assert len(df)
            df = df.rename(columns={"p": "p input to binomial"})
            df = fdr_control_benjamini_hochberg(
                df, "p-value binomial", "benjamini binomial"
            )
            df = fdr_control_benjamini_hochberg(
                df, "p-value hypergeometric", "benjamini hypergeometric"
            )
            columns = ["Group", "Set", "benjamini binomial", "benjamini hypergeometric"]
            for c in df.columns:
                if c not in columns:
                    columns.append(c)
            df = df[columns]
            return df

        job = ppg.CachedAttributeLoadingJob(
            self.cache_path
            / (hashlib.md5(self.name.encode("utf-8")).hexdigest() + "_calc"),
            self,
            "df",
            calc,
        )
        (
            job.calc.depends_on(self.regulatory_regions.load())
            .depends_on(self.query_gis.load())
            .depends_on(self.non_gap_regions.load())
            .depends_on(_add_regulatory_regions_by_gene_to_gr(self.regulatory_regions))
        )
        return job

    def calculate_regulator_region_size(self, genes):
        """Given a set of genes, calculate the size of the combined regulatory region"""
        # the problem is that the regulatory regions of genes might overlap
        # so we have to take care of that ;).
        # TODO: Replace with nested intervals.
        from mbf_nested_intervals import IntervalSet

        self.non_gap_regions.do_build_intervals()

        regions_by_chr = {}
        for g in genes:
            try:
                for (
                    chr,
                    start,
                    stop,
                ) in self.regulatory_regions.regulatory_region_by_gene[g]:
                    # TODO: continue here
                    if chr not in regions_by_chr:
                        regions_by_chr[chr] = []
                    regions_by_chr[chr].append((start, stop))
            except KeyError:
                pass  # that gene in the set is not actually in the genome...
        covered_bases = 0
        for chr, regs in regions_by_chr.items():
            non_gaps = self.non_gap_regions._interval_sets[chr]
            iv = IntervalSet.from_tuples(regs)
            iv = iv.filter_to_overlapping_and_split(non_gaps)
            covered_bases += iv.covered_units()
        return covered_bases

        for chr in regions_by_chr:
            regions_by_chr[chr].sort()
            last_stop = 0
            for start, stop in regions_by_chr[chr]:
                # this region is completly contained (sorted starts!) by the
                # one before
                if stop < last_stop:
                    continue
                if (
                    start < last_stop
                ):  # this region overlaps the one before, but ends after...
                    start = last_stop
                covered_bases += stop - start
                last_stop = stop
        return covered_bases

    def find_genes_hit_by_query(self):
        """Finds [(gene_id, interval_id)] for our query GIS"""
        hit_genes = list()
        interval_id = 0
        for dummy_idx, row in self.query_gis.df[["chr", "start", "stop"]].iterrows():
            chr, interval_start, interval_stop = row["chr"], row["start"], row["stop"]
            interval_id += 1
            if not self.non_gap_regions.has_overlapping(
                chr, interval_start, interval_stop
            ):
                raise ValueError(
                    "An interval was outside of the genomes non-gap region (%s). That makes no sense. Interval in question was %s %i %i, overlap is %s"
                    % (
                        self.non_gap_regions.name,
                        chr,
                        interval_start,
                        interval_stop,
                        self.non_gap_regions.get_overlapping(
                            chr, interval_start, interval_stop
                        ),
                    )
                )

            for dummy_idx, row in self.regulatory_regions.get_overlapping(
                chr, interval_start, interval_stop
            ).iterrows():
                assert isinstance(row["name"], str)
                hit_genes.append((row["name"], interval_id))

        return hit_genes

    def write(self, output_filename=None):
        """Write results with any hits to an output file"""
        if not output_filename:
            output_filename = Path(
                self.query_gis.result_dir, "functional annotation with GREAT.xls"
            )

        def do_write(output_filename):
            # df = self.df[self.df["Hits"] > 0].sort_values("benjamini binomial")
            df = self.df
            df = df.sort_values("benjamini binomial")
            df = df[df["benjamini binomial"] < 0.05]
            output_filename.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_filename, sep="\t")

        return ppg.FileGeneratingJob(output_filename, do_write).depends_on(
            self.perform().load
        )


# genome covered bases.
# 2858034764
# 2858034764

# all regulatory domains
# 60732793
# 60752793

# p values...
# 6.964116138947001e-09
# == 6.906703e-09
