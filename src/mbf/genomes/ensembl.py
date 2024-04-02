import re
from pathlib import Path
import pandas as pd
import mbf.externals
from mbf.externals.util import (
    # download_file_and_gzip,
    download_file,
    lazy_property,
    get_page,
    lazy_method,
)
from .base import (
    GenomeBase,
    GenomePrebuildMixin,
    include_in_downloads,
    class_with_downloads,
    MsgPackProperty,
    msgpack_unpacking_class,
)
import pypipegraph as ppg
from .common import EukaryoticCode
import pkg_resources


def download_gunzip_and_attach(url, unzipped_filename, files_to_attach):
    import shutil
    import gzip
    import tempfile

    tf = tempfile.NamedTemporaryFile(suffix=".gz")
    download_file(url, tf)
    tf.flush()

    attach = b""
    for f in files_to_attach:
        attach += Path(f).read_bytes()

    with gzip.GzipFile(tf.name, "rb") as gz_in:
        with open(unzipped_filename, "wb") as op:
            shutil.copyfileobj(gz_in, op)
            op.write(attach)


_ensembl_genome_cache = {}
data_path = Path(pkg_resources.resource_filename("mbf.genomes", "data/"))


def EnsemblGenome(species, revision, prebuild_manager=None, do_download=True):
    if not isinstance(revision, int):
        revision = int(revision)
    if prebuild_manager is None:  # pragma: no cover
        prebuild_manager = mbf.externals.get_global_manager()
    if ppg.util.global_pipegraph is not None:
        if not hasattr(ppg.util.global_pipegraph, "_ensembl_genome_dedup"):
            ppg.util.global_pipegraph._ensembl_genome_dedup = {}
        cache = ppg.util.global_pipegraph._ensembl_genome_dedup
    else:
        cache = _ensembl_genome_cache
    if (species, revision) in cache:
        res = cache[species, revision]
        if res.prebuild_manager != prebuild_manager:  # pragma: no cover
            raise ValueError(
                "Changing prebuild manager within one pipegraph is not supported"
            )
        return res
    else:
        res = _EnsemblGenome(species, revision, prebuild_manager, do_download)
        cache[species, revision] = res
        return res


@msgpack_unpacking_class
@class_with_downloads
class _EnsemblGenome(GenomePrebuildMixin, GenomeBase):
    def __init__(self, species, revision, prebuild_manager, do_download=True):
        super().__init__()
        self.prebuild_manager = prebuild_manager

        self.species = species
        if not re.match(r"^[A-Z][a-z]+_[a-z]+$", species):
            raise ValueError("Species must be capitalized like 'Homo_sapiens")
        self.revision = str(int(revision))
        self.name = f"{self.species}_{self.revision}"

        if ppg.inside_ppg():
            ppg.util.assert_uniqueness_of_object(self)

        if hasattr(
            ppg, "is_ppg2"
        ):  # SharedMultiFileGeneartingJobs will sucessfully keep them apart only if necessary
            self.prebuild_prefix = f"ensembl/{self.species}"
            if ppg.inside_ppg():
                # with ppg2, we must have only one revision per ppg
                n = self.name
                self.name = self.species
                try:
                    ppg.util.assert_uniqueness_of_object(self)
                except ValueError:
                    raise ValueError(
                        "In ppg2, you can only use one revision of a species per ppg."
                    )
                self.name = n

        else:
            self.prebuild_prefix = f"ensembl/{self.species}_{self.revision}"  # ppg1 prebuild jobs need to seperate them by version

        self.genetic_code = EukaryoticCode
        if do_download:
            self.download_genome()
        self._seq_region_is_canonical = {}
        self._canonical_cache = {}

    def __repr__(self):
        return f"EnsemblGenome({self.species}, {self.revision})"

    @include_in_downloads
    def _pb_find_server(self):
        ensembl_urls = [
            "ftp://ftp.ensembl.org/pub/release-%i/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/fungi/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/metazoa/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/plants/",
            "ftp://ftp.ensemblgenomes.org/pub/release-%i/protists/",
            # "http://ftp.ensemblgenomes.org/pub/release-%i/bacteria/", # bacteria are complicated / subdivided?
        ]

        def find_ensembl_server(output_path):
            for proto_url in ensembl_urls:
                url = proto_url % (int(self.revision),) + "fasta/"
                r = get_page(url)
                if self.species.lower() in r:
                    (output_path / "url.txt").write_text(
                        proto_url % (int(self.revision),)
                    )
                    return
            raise ValueError("Could not find this species on any ensembl server")

        server_job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/server",
            # we don't use the version for this, since we need it for building
            # various aligner versioned indices
            "1",
            [],
            ["url.txt"],
            find_ensembl_server,
            minimum_acceptable_version="1",
            maximum_acceptable_version="1",
        )
        if ppg.inside_ppg() and hasattr(ppg, "is_ppg2"):
            server_job.depends_on(
                ppg.ParameterInvariant(f"Ensembl{self.species}_revision", self.revision)
            )
        self.server_job = server_job
        return server_job

    @include_in_downloads
    def _pb_download_gtf(self):
        return self._pb_download_and_gunzip(
            "gtf",
            "gtf/" + self.species.lower() + "/",
            (rf"{self.species}\..+\.\d+\.gtf\.gz",),
            "genes.gtf",
        )  # can't have this unziped star wants it unziped

    def get_additional_gene_gtfs(self):
        if self.species == "Homo_sapiens":
            if int(self.revision) <= 74:
                return [data_path / "ribosomal_genes_grch37.gtf.gz.full.gtf.gz"]
            else:
                return [data_path / "ribosomal_genes_grch38.gtf.gz.full.gtf.gz"]
        elif self.species == "Mus_musculus":
            if int(self.revision) > 67:
                return [data_path / "ribosomal_genes_mm10.gtf.gz.full.gtf.gz"]
        elif self.species == "Ustilago_maydis":
            return [data_path / "ustilago_maydis_a2_locus.gff"]
        return []

    def get_additional_fastas(self):
        """Add additional fasta files to the genome.

        They are considered true chromosomes if 'chromosome: something'
        is in the fasta description line.
        """

        if self.species == "Ustilago_maydis":
            return [data_path / "ustilago_maydis_a2_locus.fasta"]
        return None

    @property
    def gene_gtf_dependencies(self):
        return [self._pb_download_gtf()]

    @include_in_downloads
    @lazy_method
    def _pb_download_genome_fasta(self):
        additional_fastas = self.get_additional_fastas()
        if additional_fastas:
            return self._pb_download(
                pb_name="dna",
                url="fasta/" + self.species.lower() + "/dna/",
                regexps=(
                    rf"{self.species}\..+\.dna.primary_assembly.fa.gz",
                    rf"{self.species}\..+\.dna.toplevel.fa.gz",
                ),
                output_filename="genome.fasta",
                download_func=lambda url, unzipped_filename: download_gunzip_and_attach(
                    url, unzipped_filename, additional_fastas
                ),
                additional_input_files=additional_fastas,
            )

        else:
            return self._pb_download_and_gunzip(
                "dna",
                "fasta/" + self.species.lower() + "/dna/",
                (
                    rf"{self.species}\..+\.dna.primary_assembly.fa.gz",
                    rf"{self.species}\..+\.dna.toplevel.fa.gz",
                ),
                "genome.fasta",
            )

    @include_in_downloads
    def _pb_extract_keys_from_genome(self):
        output_filename = "references.txt"

        def extract(output_path):
            from .common import iter_fasta

            fn = self.find_file("genome.fasta")
            keys = []
            for key, seq in iter_fasta(fn):
                keys.append(key)
            (output_path / output_filename).write_bytes(b"\n".join(keys))

        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/chromosomes_and_contigs",
            "1",
            [],
            [output_filename],
            extract,
        ).depends_on(self._pb_download_genome_fasta())
        return job

    @include_in_downloads
    def _pb_download_cdna_fasta(self):
        return self._pb_download_and_gunzip(
            "cdna",
            "fasta/" + self.species.lower() + "/cdna/",
            (rf"{self.species}\..+\.cdna.all.fa.gz",),
            "cdna.fasta",
        )

    @include_in_downloads
    def _pb_download_protein_fasta(self):
        return self._pb_download_and_gunzip(
            "pep",
            f"fasta/{self.species.lower()}/pep/",
            (rf"{self.species}\..+\.pep.all.fa.gz",),
            "pep.fasta",
        )

    @include_in_downloads
    def _pb_download_sql_table_definitions(self):
        return self._pb_download_straight(
            "sql/core/sql_def",
            "mysql/",
            (rf"{self.species.lower()}_core_.+",),
            "core.sql.gz",
            lambda match: f"{match.strip()}/{match.strip()}.sql.gz",
        )

    def _pb_download_sql_table(self, table_name):
        """Helper to download sql tables as mysql dumps"""
        job = self._pb_download_straight(
            f"sql/core/{table_name}",
            "mysql/",
            (rf"{self.species.lower()}_core_.+",),
            f"{table_name}.txt.gz",
            lambda match: f"{match.strip()}/{table_name.strip()}.txt.gz",
        ).depends_on(self._pb_download_sql_table_definitions())
        job.table_name = table_name
        return job

    @include_in_downloads
    @lazy_method
    def _pb_download_sql_tables(self):
        tables = [
            ("gene"),  # for description
            ("transcript"),  # for external name lookup transcript -> gene
            ("translation"),  # for external name lookup translation -> gene
            ("stable_id_event"),  # for stable_id changes
            ("external_db"),  # for external name lookup
            ("object_xref"),  # for external name lookup
            ("xref"),  # for external name lookup
            ("alt_allele"),  # for finding 'canonical' ids
            ("seq_region"),  # for finding 'canonical' ids
            ("seq_region_attrib"),  # for finding 'canonical' ids
            ("attrib_type"),  # for finding 'canonical' ids
            ("external_synonym"),  # for finding gene names...
        ]
        return [self._pb_download_sql_table(x) for x in tables]

    # def _pb_download_and_gzip(self, pb_name, url, regexps, output_filename):
    # return self._pb_download(
    # pb_name, url, regexps, output_filename, download_file_and_gzip
    # )

    @lazy_property
    def base_url(self):
        return self.server_job.find_file("url.txt").read_text()

    @lazy_method
    def get_true_chromosomes(self):
        """Get the names of 'true' chromosomes, ie. no scaffolds/contigs
        in genomes that have chromosomes, otherwise all"""
        fn = self.find_file("references.txt")
        keys = Path(fn).read_text().split("\n")
        chroms = [x for x in keys if "chromosome:" in x]
        if not chroms:
            chroms = keys
        return [x[: x.find(" ")] for x in chroms]

    def _load_from_sql(
        self, table_name, columns=None, check_for_columns=None, **kwargs
    ):
        table_columns = self._get_sql_table_column_names(table_name)
        for c in columns:
            if not c in table_columns:
                raise ValueError(c, "available", table_columns)

        df = pd.read_csv(
            self.find_file(f"{table_name}.txt.gz"),
            sep="\t",
            header=None,
            names=table_columns,
            usecols=columns,
            na_values="\\N",
            lineterminator="\n",
            escapechar="\\",
            **kwargs,
        )
        if check_for_columns:
            for c in check_for_columns:
                if not c in table_columns:  # pragma: no cover
                    raise KeyError(c, "availabel", table_columns)
        return df

    def _prepare_df_genes_meta(self):
        """Meta data for genes.
        Currently contains:
            'description'
        """
        try:
            df = self._load_from_sql(
                "gene", ["stable_id", "description"], ["stable_id"]
            )
        except KeyError:  # pragma: no cover
            raise ValueError(
                "No stable_id column found - "
                "old ensembl, split into seperate table, add support code?"
            )
        res = df.set_index("stable_id")
        res.index.name = "gene_stable_id"
        return res

    df_genes_meta = MsgPackProperty(
        lambda self: [
            x for x in self._pb_download_sql_tables() if x.table_name == "gene"
        ]
    )

    def _get_sql_table_column_names(self, sql_table_name):
        """Read the sql definition and extract column names"""
        import gzip

        with gzip.GzipFile(self.find_file("core.sql.gz")) as op:
            raw = op.read().decode("utf-8")
            if "\n-- Table structure" in raw:
                parts = raw.split("\n-- Table structure")[1:]
            else:
                parts = raw.split("\n\n")
            for p in parts:
                if "CREATE TABLE" not in p:
                    continue
                p = p[p.find("CREATE TABLE") :]
                if "  PRIMARY" in p:
                    p = p[: p.find("  PRIMARY")]
                elif "  UNIQUE" in p:
                    p = p[: p.find("  UNIQUE")]
                elif "  KEY" in p:
                    p = p[: p.find("  KEY")]

                else:  # pragma: no cover
                    raise ValueError(p)
                names = re.findall("`([^`]+)`", p)
                table_name, *columns = names
                if table_name == sql_table_name:
                    return columns
        raise KeyError(f"{sql_table_name} not in core.sql.gz")  # pragma: no cover

    def _prepare_lookup_stable_id_events(self):
        """Lookup old_stable_id -> new_stable_id"""
        df = self._load_from_sql("stable_id_event", ["old_stable_id", "new_stable_id"])
        lookup = {}
        olds = [str(x) for x in df["old_stable_id"].values]
        news = [str(x) for x in df["new_stable_id"].values]
        for old in olds:
            lookup[old] = set()
        for old, new in zip(olds, news):
            lookup[old].add(new)
        return pd.DataFrame(
            {"old": list(lookup.keys()), "new": [list(x) for x in lookup.values()]}
        ).set_index("old")

    lookup_stable_id_events = MsgPackProperty(
        lambda self: [
            x
            for x in self._pb_download_sql_tables()
            if x.table_name == "stable_id_event"
        ]
    )

    def newest_stable_ids_for(self, stable_id):
        """Get the most up to date and current stable_ids for genes, transcripts, proteins).
        Plurarl for gene might have split, or have been deleted.
        returns a set of new ids.
        """
        try:
            valid_ids = set(self.df_genes.index)
            valid_ids.update(self.df_transcripts.index)
            valid_ids.update(self.df_proteins.index)
            res = set(self.lookup_stable_id_events.loc[stable_id]["new"])
            res = set(
                [x for x in res if x in valid_ids]
            )  # filter those that are no longer in the database - no matter that they were m             apped somewhere else in between
            return res
        except KeyError as e:
            # see if it's a current id where we're simply lacking the stable_id_event for some reason
            if stable_id in valid_ids:
                return set([stable_id])
            else:
                raise e

    def get_external_dbs(self):
        """Return the names of all external dbs that actually have xrefs"""
        df_external_db = self._load_from_sql(
            "external_db", ["external_db_id", "db_name"]
        )
        with_data = set(
            self._load_from_sql("xref", ["external_db_id"])["external_db_id"].unique()
        )
        return sorted(
            df_external_db["db_name"][df_external_db["external_db_id"].isin(with_data)]
        )

    def get_external_db_to_gene_id_mapping(self, external_db_name):
        """Return a dict external id -> set(stable_id, ...)
        for a given external db - e.g. EntrezGene, list
        with get_external_dbs()
        """
        df_external_db = self._load_from_sql(
            "external_db", ["external_db_id", "db_name"]
        ).set_index("db_name")
        external_db_id = df_external_db.at[external_db_name, "external_db_id"]
        xref = self._load_from_sql(
            "xref", ["dbprimary_acc", "external_db_id", "xref_id"]
        ).set_index("xref_id")
        xref = xref[xref.external_db_id == external_db_id]
        object_xref = self._load_from_sql(
            "object_xref", ["ensembl_object_type", "xref_id", "ensembl_id"]
        )
        object_xref = object_xref[object_xref.xref_id.isin(set(xref.index))]
        # object_xref = object_xref[object_xref["ensembl_object_type"] == "Gene"]
        result = {}
        transcripts = None
        translations = None
        genes = self._load_from_sql("gene", ["gene_id", "stable_id"]).set_index(
            "gene_id"
        )
        for row in object_xref.itertuples(index=False):
            if row.ensembl_object_type == "Gene":
                gene_id = row.ensembl_id
            elif row.ensembl_object_type == "Transcript":
                if transcripts is None:
                    transcripts = self._load_from_sql(
                        "transcript", ["transcript_id", "gene_id"]
                    ).set_index("transcript_id")

                gene_id = transcripts.at[row.ensembl_id, "gene_id"]
            elif row.ensembl_object_type == "Translation":
                if translations is None:
                    translations = self._load_from_sql(
                        "translation", ["translation_id", "transcript_id"]
                    ).set_index("translation_id")
                if transcripts is None:
                    transcripts = self._load_from_sql(
                        "transcript", ["transcript_id", "gene_id"]
                    ).set_index("transcript_id")

                transcript_id = translations.at[row.ensembl_id, "transcript_id"]
                gene_id = transcripts.at[transcript_id, "gene_id"]
            else:
                print(row)
                raise ValueError("Mapped to neiter a transcript, nor a gene")

            gene_stable_id = genes.at[gene_id, "stable_id"]

            db_primary = xref.at[row.xref_id, "dbprimary_acc"]
            if not db_primary in result:
                result[db_primary] = set()
            result[db_primary].add(gene_stable_id)
        return result

    def get_external_db_to_translation_id_mapping(self, external_db_name):
        """Return a dict external id -> set(translation_stable_id, ...)
        for a given external db - e.g. Uniprot/SWISSPROT
        see get_external_dbs() for a list
        """
        df_external_db = self._load_from_sql(
            "external_db", ["external_db_id", "db_name"]
        ).set_index("db_name")
        external_db_id = df_external_db.at[external_db_name, "external_db_id"]
        xref = self._load_from_sql(
            "xref", ["dbprimary_acc", "external_db_id", "xref_id"]
        ).set_index("xref_id")
        xref = xref[xref.external_db_id == external_db_id]
        object_xref = self._load_from_sql(
            "object_xref", ["ensembl_object_type", "xref_id", "ensembl_id"]
        )
        object_xref = object_xref[object_xref.xref_id.isin(set(xref.index))]
        # object_xref = object_xref[object_xref["ensembl_object_type"] == "Gene"]
        result = {}
        translations = self._load_from_sql(
            "translation", ["translation_id", "stable_id", "version"]
        ).set_index("translation_id")

        for row in object_xref.itertuples(index=False):
            if row.ensembl_object_type == "Translation":
                translation_stable_id = (
                    translations.at[row.ensembl_id, "stable_id"]
                    + "."
                    + str(translations.at[row.ensembl_id, "version"])
                )
            else:
                print(row)
                raise ValueError("not at translation mapping")

            db_primary = xref.at[row.xref_id, "dbprimary_acc"]
            if not db_primary in result:
                result[db_primary] = set()
            result[db_primary].add(translation_stable_id)
        return result

    @lazy_property
    def allele_groups(self):
        df = self._load_from_sql("alt_allele", ["alt_allele_group_id", "gene_id"])
        gene_df = self._load_from_sql("gene", ["gene_id", "stable_id"]).rename(
            columns={"gene_stable_id": "stable_id"}
        )
        df = df.join(gene_df.set_index("gene_id"), "gene_id")
        return df.set_index("stable_id")

    def name_to_canonical_id(self, name, break_ties_by_number_of_transcripts=False):
        """Given a gene name, lookup up it's stable ids, and return the
        one that's on the primary assembly from the allele group"""
        key = name, break_ties_by_number_of_transcripts
        if not key in self._canonical_cache:
            r = self._name_to_canonical_id(name, break_ties_by_number_of_transcripts)
            self._canonical_cache[key] = r
        else:
            r = self._canonical_cache[name, break_ties_by_number_of_transcripts]
        return r

    def _name_to_canonical_id(self, name, break_ties_by_number_of_transcripts=False):
        name_candidates = set(
            [x for x in self.name_to_gene_ids(name) if not x.startswith("LRG")]
        )
        if not name_candidates:  # pragma: no cover
            raise KeyError("No gene named %s" % name)
        ag = self.allele_groups
        ag_ids = [
            x
            for x in ag.alt_allele_group_id[
                ag.alt_allele_group_id.index.isin(name_candidates)
            ].unique()
            if not pd.isnull(x)
        ]
        ag_candidates = set(ag.index[ag.alt_allele_group_id.isin(ag_ids)])
        if len(ag_ids) == 1 and name_candidates.issubset(ag_candidates):
            # the easy case, everything matches
            on_primary = [
                x
                for x in ag_candidates
                if x
                in self.df_genes.index  # for there is no entry in genes.gtf if it's not on a not 'non_ref' chromosome.
            ]
            if len(on_primary) == 1:
                return on_primary[0]
            elif len(on_primary) == 0:  # pragma: no cover
                # if self.species == "Homo_sapiens" and name == "HLA-DRB3": # HLA-DRB3 is not in genes.gtf!
                # known issue - return basically any of the candidates on alternate regions, but be consistent.
                # return sorted(ag_candidates)[0]
                raise ValueError("No primary gene found for %s" % name)
            else:  # pragma: no cover
                raise ValueError(
                    "Multiple gene on primary assemblies found for %s" % name
                )
        elif len(ag_ids) == 0 and len(name_candidates) == 1:  # pragma: no cover
            # another easy case, there are no alternatives
            return list(name_candidates)[0]
        else:
            if break_ties_by_number_of_transcripts:
                name_candidates = list(name_candidates)
                name_candidates.sort(
                    key=lambda gene_stable_id: len(
                        self.genes[gene_stable_id].transcripts
                    )
                )
                return name_candidates[-1]
            else:
                raise ValueError(  # pragma: no cover
                    "Could not determine canonical gene for '%s'. "
                    "Either pass break_ties_by_number_of_transcripts=True, "
                    "or use name_to_gene_ids()"
                    " and have a look yourself (don't forget the allele groups).\n"
                    "Name candidates: %s\n"
                    "AG candidates: %s\n"
                    "AG ids: %s" % (name, name_candidates, ag_candidates, ag_ids)
                )

    def load_synonyms(self):
        cn = "_load_synonyms"
        if not hasattr(self, cn):
            try:
                synonyms = self._load_from_sql(
                    "external_synonym", ["xref_id", "synonym"]
                ).set_index("xref_id")
                object_xref = self._load_from_sql(
                    "object_xref", ["ensembl_id", "ensembl_object_type", "xref_id"]
                )
                object_xref = (
                    object_xref[object_xref.ensembl_object_type == "Gene"]
                    .set_index("ensembl_id")
                    .drop(columns="ensembl_object_type")
                )
                genes = self._load_from_sql("gene", ["stable_id", "gene_id"])
                og = genes.join(object_xref, "gene_id").set_index("xref_id")
                s = synonyms.join(og, "xref_id")[["synonym", "stable_id"]]
                s = s.assign(synonym=s.synonym.str.lower())
                setattr(self, cn, s.set_index("synonym"))
            except Exception as e:
                raise ValueError(e)
        return getattr(self, cn)

    def synonym_to_gene_ids(self, synonym):
        return list(self.load_synonyms().loc[[synonym.lower()]].stable_id)

    def load_orthologes(self, target_species):
        if not hasattr(self, "_orthologes"):
            self._orthologes = {}
        if not target_species in self._orthologes:
            from .ensembl_compara import EnsemblCompara

            compara = EnsemblCompara(int(self.revision))
            input_df = compara.get_homology_table(self.species, target_species)
            # print(input_df.columns)
            self._orthologes[target_species] = {}
            for homology_type, sub_df in (
                input_df[
                    [
                        "gene_stable_id",
                        "homology_type",
                        "homology_gene_stable_id",
                    ]
                ]
                .sort_values("homology_type")
                .groupby("homology_type")
            ):
                r = {}
                for from_, to_ in zip(
                    sub_df["gene_stable_id"], sub_df["homology_gene_stable_id"]
                ):
                    if not from_ in r:
                        r[from_] = set()
                    r[from_].add(to_)
                self._orthologes[target_species][homology_type] = r
        return self._orthologes[target_species]

    def to_orthologs(
        self,
        gene_stable_id,
        target_species,
        types=["ortholog_one2one", "ortholog_one2many", "ortholog_many2many"],
    ):
        """Convert a gene stable id to it's orthologues"""
        self.load_orthologes(target_species)
        res = set()
        for t in types:
            res.update(self._orthologes[target_species][t].get(gene_stable_id, []))
        return res

    @property
    def assembly(self):
        if self.species == "Homo_sapiens" and int(self.revision) >= 96:
            return "hg38"
        elif self.species == "Mus_musculus" and int(self.revision) > 67:  # I believe...
            return "mm10"
        else:
            raise ValueError("fNo assembly name for {self}, update the function?")

