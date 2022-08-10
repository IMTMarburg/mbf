"""This module provides support for generating your own genesets as well as
interfaces to serveral databases of such.

A GeneGroup is an object that supplies dependencies, and can return genesets
that have been tailored to a genome.  (genesets are dictionaries of {'name':
set([stable_id1, stable_id2...,])}).

"""

import os
import pypipegraph as ppg

import mbf_genomes
from kitchen.text.converters import to_unicode
from ordered_set import OrderedSet
import re
from mbf_fileformats.util import open_file
from pathlib import Path
from mbf_externals.util import lazy_method


class _GroupsBase(object):
    """This class parses the following file format:
    Fileformat:
    # This is a comment line. It is ignored everywhere
    !parse_mode Name of Geneset
    #that's a space or a tab after parse_mode
    #commas and spaces are treated interchangably
    gene1 gene2,gene3
    gene4
    gene5
    gene6 gene7
    #and on to the next gene set.
    !parse_mode name of next geneset
    gene3
    gene4
    gene8

    Valid parse_modes are
        stable_id - ensembl stable ids of any species
            HUGO - http://www.genenames.org/
            mouse_name - mouse gene names


        Genes are automagically translated via ensembl compara if the species
        from the file and the one requested do not match.
    """

    def __init__(self, group_name, path):
        self.path = path
        self.name = group_name
        self.cached_gene_stable_ids = None
        self.cached_gene_stable_ids_for_genome = None

    def get_sets(self, genome):
        """Return the sets known for this group in the given genome"""
        raise NotImplementedError()

    def get_dependencies(self):
        """Do the sets dependen on anything (pypipegraph)?"""
        raise NotImplementedError()

    def get_url(self, name):
        """Groups could return an url for their sets here that offers further information"""
        return ""

    def parse(self, filename, genome):
        """Parse the file into an iterator yielding (name, gene_stable_ids) (in the target genome."""
        handle = open(filename, "r")
        mode = False
        name = False

        def translate_genes():
            """depending on the mode, dispatch to the right function"""
            genes = False
            if mode == "HUGO":  # HUGO http://www.genenames.org/
                was_ok, genes_human = self.translate_hugo_human(
                    genes_in_file_format, genome
                )
                if not was_ok:
                    raise ValueError(
                        "Could not translate HUGOs to stable_id - dataset %s:\n %s"
                        % (filename, "\n".join(genes_human))
                    )
                if genome.species != "Homo_sapiens":
                    genes = self.convert_species("Homo_sapiens", genome, genes_human)
                else:
                    genes = genes_human
            elif mode == "stable_id":
                was_ok, genes = self.translate_stable_ids(genes_in_file_format, genome)
                if not was_ok:
                    raise ValueError(
                        "Could not map the following stable_ids - dataset %s:\n %s"
                        % (filename, "\n".join(genes))
                    )
            elif mode == "mouse_name":  # mouse gene names
                was_ok, genes_mouse = self.translate_mouse_gene_names(
                    genes_in_file_format, genome
                )
                if not was_ok:
                    raise ValueError(
                        "Could not map the following mouse gene names - dataset %s:\n %s"
                        % (filename, "\n".join(genes_mouse))
                    )
                if genome.species != "Mus_musculus":
                    genes = self.convert_species("Homo_sapiens", genome, genes_human)
                else:
                    genes = genes_mouse
            elif mode:
                raise ValueError("Could not interpret mode: %s" % mode)
            if genes:
                genes = OrderedSet(genes)
            return genes  # we can return false if mode was not set yet

        for line_no, row in enumerate(handle.readlines()):
            row = row.strip()
            if not row or row[0] == "#":  # ignore comment lines or empty lines
                continue
            elif row.startswith("!"):
                genes = translate_genes()
                if genes:
                    yield name, genes
                genes_in_file_format = []
                mode = row.split()[0][1:]  # take space and tab.
                name = " ".join(row.split()[1:])
                if not name:
                    raise ValueError(
                        "No name set in line %i of %s" % (line_no, filename)
                    )
                name = to_unicode(name, "utf-8")
            elif not mode:
                raise ValueError(
                    "File did not contain a !mode name line before genes: Line: %i %s"
                    % (line_no, filename)
                )
            elif mode:
                row = row.replace(",", " ")
                genes_in_file_format.extend(row.split())
            else:
                raise ValueError("Should not be reached")
        genes = translate_genes()
        if genes:
            yield name, genes
        handle.close()

    def translate_hugo_human(self, hugo_ids, genome, ignore_errors=False):
        """HUGO gene symbols - see  http://www.genenames.org/"""
        hs_rev = int(genome.revision)
        human_genome = mbf_genomes.EnsemblGenome("Homo_sapiens", hs_rev)

        if not hasattr(
            self, "cached_gene_stable_ids_for_genome"
        ) or self.cached_gene_stable_ids_for_genome != (
            "Homo_sapiens",
            genome.revision,
        ):
            self.cached_gene_stable_ids_for_genome = ("Homo_sapiens", genome.revision)
            self.cached_gene_stable_ids = set(human_genome.df_genes.index.values)
        gene_stable_ids = self.cached_gene_stable_ids
        # lookup = human_genome.get_hugo_to_stable_id_lookup()
        genes = []
        errors = []
        for hugo_name in hugo_ids:
            stable_ids = False
            try:
                # stable_ids = list(lookup[hugo_name.upper()]) #Hugo is always upper
                if hugo_name in gene_stable_ids:
                    stable_ids = [hugo_name]
                else:
                    stable_ids = human_genome.name_to_gene_ids(hugo_name)
                for x in stable_ids:
                    if "," in x:
                        raise ValueError(
                            "hugo name lookup returned a comma separated list. Reintroduce the splitting suggested by marco"
                        )
                        # stable_ids = re.split(',\\s+', lookup[hugo_name.upper()]) #Hugo is always upper
            except KeyError:
                stable_ids = list(human_genome.alternative_name_to_gene_ids(hugo_name))
            if stable_ids:
                genes.extend(stable_ids)
        if errors and not ignore_errors:
            return False, errors
        else:
            return True, genes

    def translate_mouse_gene_names(self, hugo_ids, genome):
        """Mouse gene names, as defined by ensembl"""
        mouse_genome = mbf_genomes.EnsemblGenome("Mus_musculus", genome.revision)
        genes = []
        errors = []
        for hugo_name in hugo_ids:
            stable_ids = False
            try:
                stable_ids = list(mouse_genome.name_to_gene_ids(hugo_name))
            except KeyError:
                alternative = list(mouse_genome.alternative_name_to_gene_ids(hugo_name))
                if alternative:
                    stable_ids = alternative
                else:
                    errors.append(hugo_name)
            if stable_ids:
                genes.extend(stable_ids)
        if errors:
            return False, errors
        return True, genes

    def translate_stable_ids(self, stable_ids, target_genome, check_errors=True):
        """Straight ensembl stable_ids of any species"""
        first = list(stable_ids)[0]
        if first.startswith("ENSG"):
            source_species_name = "Homo_sapiens"
        elif first.startswith("ENSMUSG"):
            source_species_name = "Mus_musculus"
        else:
            raise ValueError(
                "You need to extend translate_stable_ids to discern the species from genes like %s"
                % first
            )
        source_genome = mbf_genomes.EnsemblGenome(
            source_species_name, target_genome.revision
        )
        print(source_genome)
        errors = []
        genes = []
        for s in stable_ids:
            try:
                source_genome.get_gene_info(s)  # check if it's in the genome.
                genes.append(s)
            except KeyError:
                # errors.append(s)
                try:
                    for s_new in source_genome.newest_stable_ids_for(s):
                        genes.append(s_new)
                except KeyError:
                    errors.append(s)
        if check_errors and errors:
            return False, errors
        if source_genome.species != target_genome.species:
            genes = self.convert_species(source_species_name, target_genome, genes)
        return True, genes

    def convert_species(
        self, source_species_name, target_genome, gene_stable_ids,
    ):
        """Use Ensembl Compara for homology deduction"""
        source_genome = mbf_genomes.EnsemblGenome(
            source_species_name, int(target_genome.revision)
        )

        converted_genes = set()
        for source_stable_id in gene_stable_ids:
            homologues = source_genome.to_orthologs(
                source_stable_id, target_genome.species
            )
            converted_genes.update(
                homologues
            )  # we keep all possible mappings. The right thing to do for functional groups, I'd argue
        return list(converted_genes)

    def _set_name_to_group_id(self, set_name):
        return set_name


class GroupsFromDirectory(_GroupsBase):
    """Read a whole directory of files"""

    def get_sets(self, genome):
        res = {}
        for filename in os.listdir(self.path):
            for name, genes in self.parse(os.path.join(self.path, filename), genome):
                res[name] = genes
        return res

    def get_dependencies(self):
        res = []
        for filename in os.listdir(self.path):
            if filename.endswith(".txt"):
                res.append(ppg.FileTimeInvariant(os.path.join(self.path, filename)))
        return res


class GroupsFromFile(_GroupsBase):
    """Read a single file"""

    def get_sets(self, genome):
        res = {}
        for name, genes in self.parse(self.path, genome):
            res[name] = genes
        return res

    def get_dependencies(self):
        res = []
        res.append(ppg.FileTimeInvariant(self.path))
        return res


class GroupsFromFlybaseGO(_GroupsBase):
    """Read a flybase.org GO file after download from http://flybase.org/static_pages/downloads/current/go/"""

    def parse(self, filename, genome):
        import pandas

        df = pandas.read_csv(
            filename,
            sep="\t",
            skiprows=5,
            header=None,
            names=[
                "A",
                "stable_id",
                "name",
                "empty",
                "GO_term",
                "F",
                "G",
                "H",
                "I",
                "J",
                "K",
                "biotype",
                "taxon",
                "N",
                "source",
                "P",
                "Q",
            ],
        )
        # df.columns = ['A', 'stable_id', 'name', 'empty', 'GO_term', 'F', 'G', 'H', 'I', 'J', 'K','biotype', 'taxon', 'N', 'source', 'P', 'Q']
        for x in df.groupby("GO_term"):
            yield (x[0], x[1]["stable_id"].values)

    def get_sets(self, genome):
        res = {}
        for filename in os.listdir(self.path):
            if filename.endswith(".fb"):
                for name, genes in self.parse(
                    os.path.join(self.path, filename), genome
                ):
                    res[name] = genes
        return res

    def get_dependencies(self):
        res = []
        for filename in os.listdir(self.path):
            if filename.endswith(".fb"):
                res.append(ppg.FileTimeInvariant(os.path.join(self.path, filename)))
        return res


class GroupsFromGenes(_GroupsBase):
    """Read geneset from genes object"""

    def __init__(self, group_name, list_of_genes, mode="stable_id"):
        self.name = group_name
        self.mode = mode
        self.list_of_genes = list_of_genes

    def get_sets(self, genome):
        res = {}
        for name, genes in self.parse(genome):
            res[name] = genes
        return res

    def get_dependencies(self):
        res = []
        for genes in self.list_of_genes:
            res.append(genes.load())
        return res

    def parse(self, genome):
        """Parse the file into an iterator yielding (name, gene_stable_ids) (in the target genome."""
        if self.mode == "HUGO":
            column_name = "name"
        elif self.mode == "stable_id":
            column_name = "stable_id"
        for genes in self.list_of_genes:
            yield genes.name, genes.df.gcv(column_name)

    def _set_name_to_group_id(self, set_name):
        return set_name


class GenomicRegionDerivedGeneSet_ClosestTwoGenes:
    """A handy adapter to define additional reference sets for HypergeometricEnrichment
    (ie. these go into the groups, they don't work as a query)

    """

    def __init__(self, source):
        self.source = source
        self.name = "GRDGS_ClosestTwoGenes_" + source.name
        self.cache_path = os.path.join("cache", "GenomicRegions" + self.name)
        Path(self.cache_path).mkdir(exist_ok=True)

    def get_dependencies(self, genome):
        if genome != self.source.genome:
            raise ValueError(
                "This GenomicRegionDerivedGeneSet_ClosestTwoGenes has a different genome then is being queried"
            )
        return [self.load_data()]

    def get_url(self, name):
        return ""

    def load_data(self):
        def calc():
            result = {}
            for gr in self.source:
                genes = set()
                for row in gr.df.iter_rows():
                    next_genes = self.source.genome.get_genes_with_closest_start(
                        row["chr"], (row["start"] + row["stop"]) / 2.0, 2
                    )
                    for found_gene in next_genes:
                        genes.add(found_gene[0])  # stable_id
                result[gr.name] = genes
            return result

        def inject():
            return [gr.annotate() for gr in self.source]

        calc_job = ppg.CachedAttributeLoadingJob(
            os.path.join(self.cache_path, "calc_sets"), self, "_sets", calc
        )
        dep_injection_job = ppg.DependencyInjectionJob(
            os.path.join(self.cache_path, "inject_calc_dependencies"), inject
        ).depends_on(self.source.get_dependencies())
        calc_job.depends_on(dep_injection_job)
        return calc_job

    def get_sets(self, genome):
        if genome != self.source.genome:
            raise ValueError(
                "This GenomicRegionDerivedGeneSet_ClosestTwoGenes has a different genome then is being queried"
            )
        return self._sets


class MotifSourceDerivedGeneSets:
    """For each gene in the genome, get the 'promotor region', check whether each motif is present,
    and build a dict {motif_name: set(genes_with_that_motif_in_their_promotor)}
    """

    def __init__(self, motif_source, upstream_distance=-5000, downstream_distance=1000):
        if upstream_distance > 0:
            raise ValueError("upstream_distance must be < 0")
        if downstream_distance < 0:
            raise ValueError("downstream_distance must be > 0")
        self.source = motif_source
        self.upstream_distance = upstream_distance
        self.downstream_distance = downstream_distance
        self.name = "MSDGS_%i_%i_%s" % (
            upstream_distance,
            downstream_distance,
            self.source.name,
        )
        self.cache_path = os.path.join("cache", "functional", "MSDGS", self.name)
        exptools.ensure_path(self.cache_path)

    def get_dependencies(self, genome):
        if genome != self.source.genome:
            raise ValueError(
                "This MotifSourceDerivedGeneSets has a different genome then is being queried"
            )
        return [self.load_data()]

    def get_url(self, name):
        return ""

    def load_data(self):
        def calc():
            sets = {}
            genome = self.source.genome
            motifs = list(self.source)
            for motif in motifs:
                sets[motif.name] = set()
            # this is dog slow - needs multicore
            for transcript_info in genome.get_all_transcripts().iter_rows():
                gene_stable_id = transcript_info["gene_stable_id"]
                if transcript_info["strand"] == 1:
                    tss = transcript_info["start"]
                    start = tss + self.upstream_distance
                    stop = tss + self.downstream_distance
                else:
                    tss = transcript_info["stop"]
                    start = tss - self.downstream_distance
                    stop = tss - self.upstream_distance  # - - => +
                try:
                    seq = genome.get_sequence(transcript_info["chr"], start, stop)
                except IOError:
                    continue
                for motif in motifs:
                    threshold = motif.max_score * 0.5
                    cum_score, max_score = motif.scan(seq, threshold, False)
                    if max_score >= threshold:
                        sets[motif.name].add(gene_stable_id)
            return sets

        def inject():
            return [gr.annotate() for gr in self.source]

        calc_job = ppg.CachedAttributeLoadingJob(
            os.path.join(self.cache_path, "calc_sets"), self, "_sets", calc
        ).depends_on(self.source.get_dependencies())
        # jdep_injection_job = ppg.DependencyInjectionJob(os.path.join(self.cache_path, 'inject_calc_dependencies'), inject).depends_on(self.source.get_dependencies())
        # jcalc_job.depends_on(dep_injection_job)
        return calc_job

    def get_sets(self, genome):
        if genome != self.source.genome:
            raise ValueError(
                "This MotifSourceDerivedGeneSets has a different genome then is being queried"
            )
        return self._sets


db_dir = os.path.join(os.path.dirname(__file__), "datasets")


class GMTDataset(_GroupsBase):
    """interface for GSEA GMT dataset files (with HUGO gene symbols)"""

    def __init__(self, name, filename, human_genome=None):
        self.name = name
        self.filename = filename
        self.human_genome = human_genome
        self._cache = {}

    def get_dependencies(self):
        return ppg.FileTimeInvariant(self.filename)

    def get_sets(self, genome):
        if not genome in self._cache:
            res = {}
            for name, genes in self.parse(self.filename, genome):
                res[name] = genes
                self._cache[genome] = res
        return self._cache[genome]

    def parse(self, filename, genome):
        if self.human_genome is None:
            human_genome = genome
        else:
            human_genome = self.human_genome
        handle = open(filename, "r")
        for row in handle.readlines():
            row = row.strip()
            if not row:
                continue
            row = row.split("\t")
            name = row[0]
            name = to_unicode(name, "utf-8")
            # url = row[1]
            gene_symbols = row[2:]
            ok, stable_ids = self.translate_hugo_human(
                gene_symbols, human_genome, ignore_errors=True
            )
            if not ok:
                raise ValueError(
                    "Could not translate into stable ids: %s (from file %s, geneset %s)"
                    % (stable_ids, filename, name)
                )
            if genome.species != "Homo_sapiens":

                stable_ids = self.convert_species(
                    "Homo_sapiens", genome, stable_ids
                )
            yield name, set(stable_ids)
        handle.close()


class MSigDataset(GMTDataset):
    def __init__(self, cX, version="v7.2", subgroup=None, human_genome=None):
        """Create a Groups from the molecular signatura dabase (BROAD)
        cX is one of c1, c2, c3, c4, c5, c6
        positional, curated, motif, computational GO, oncogenic signatures
        (see http://www.broadinstitute.org/gsea/msigdb/collections.jsp)
        """
        self.name = "MSIG_%s" % cX
        if subgroup:
            self.name = "MSIG_%s_%s" % (cX, subgroup)
            self.filename = os.path.join(
                db_dir, "msig", "%s.%s.%s.symbols.gmt" % (cX, subgroup, version)
            )
        else:
            self.filename = os.path.join(
                db_dir, "msig", "%s.all.%s.symbols.gmt" % (cX, version)
            )
        self._cache = {}
        self.human_genome = human_genome

    def get_url(self, set_name):
        return "http://www.broadinstitute.org/gsea/msigdb/cards/%s" % set_name


def _parse_panther_hierarchie(filehandle):
    """Parse hierarchie from a the second column of a file looking like
    BP00001 2.01.00.00.00   Carbohydrate metabolism Metabolic processes of carbohydrates, including the breakdown and biosynthesis of carbohydrates.    GO:0005975  carbohydrate metabolism
    """
    by_key = {}
    desc_by_key = {}
    id_to_desc = {}
    for line in filehandle:
        line = line.split("\t")
        pc_id = line[0]
        key = line[1]
        by_key[key] = pc_id
        desc_by_key[key] = line[-1].strip()
        id_to_desc[pc_id] = line[2].strip()
    parents = {}
    children = {}
    for key, pc_id in by_key.items():
        id_parts = key.split(".")
        try:
            last_value = id_parts.index("00") - 1
        except ValueError:
            last_value = len(id_parts) - 1
        if last_value > 0:  # otherwise this is a top level and doesn't have parents.
            id_parts[last_value] = "00"
            parent_key = ".".join(id_parts)
            if parent_key in by_key:
                parents[pc_id] = by_key[parent_key]  # + " ! " + desc_by_key[parent_key]
                if not by_key[parent_key] in children:
                    children[by_key[parent_key]] = set()
                children[by_key[parent_key]].add(by_key[key])
            else:
                parents[pc_id] = "top_level"
        else:
            parents[pc_id] = "top_level"
    return parents, children, id_to_desc


class PantherDBDataset:
    # __metaclass__ = common.Singleton

    def __init__(self):
        self.name = "Panther"
        self.pantherdir = "/gb/imt/datasets/pantherdb_dataset"
        self.data_path = os.path.join(db_dir, "pantherdb", "7.0")
        self.structure = (
            {}
        )  # convention is (parent, children, id_to_desc), both dicts of id -> other id, respectively id -> ids. There should be one 'parent' labled 'top_level'

    def get_sets(self, genome):
        if not self.structure:
            self.parse_structure()
        if genome.species == "Mus_musculus":
            return self.get_sets_mouse(genome)
        elif genome.species == "Homo_sapiens":
            return self.get_sets_human(genome)
        else:
            raise ValueError(
                "Currently only have panther data for mouse though there is more available at the website"
            )

    def get_sets_mouse(self, genome):
        data_file = os.path.join(
            self.data_path, "PTHR7.0_Mouse.bz2"
        )  # guess you gotta adjust this if you change species
        return self._get_sets(genome, data_file, "mgi")

    def _get_sets(self, genome, data_file, mode="mgi"):
        sets = {}
        # read in the definitions...
        with open_file(data_file, "r") as op:
            for line in op:
                line = line.strip()
                if line:
                    values = line.split("\t")
                    ids = values[0].split("|")
                    if mode == "mgi":
                        if not "MGI" in values[0]:
                            raise ValueError(
                                "Line without MGI:r%s" % line
                            )  # sanity check
                        mgi = [x for x in ids if x.startswith("MGI")][
                            0
                        ]  # should only ever be one...
                        mgi = mgi[4:]  # cut off duplicate 'MGI:'
                        ensembl_stable_id = genome.alternative_name_to_gene_ids(mgi)
                    else:
                        mgi = [x for x in ids if x.startswith("ENSEMBL:")][
                            0
                        ]  # should only ever be one...
                        ensembl_stable_id = [mgi[mgi.find(":") + 1 :]]
                    for (
                        stable_id
                    ) in (
                        ensembl_stable_id
                    ):  # a clever way to say 'if you found one', since from the mgi data there really should be only one
                        self.parse_terms(values[6:], sets, stable_id)
        # sets is a dict(name -> set)
        # and we need to transport terms upward from lower groups to their parents.
        for group_name in sets:
            parent = group_name
            while True:
                try:
                    parent = self.get_parent(parent)
                    if parent and parent in sets:
                        sets.parent.update(sets[group_name])
                    else:
                        break
                except KeyError:
                    break
        return sets

    def get_sets_human(self, genome):
        data_file = os.path.join(self.data_path, "PTHR7.0_Human.bz2")
        return self._get_sets(genome, data_file, "ensembl")

    def parse_terms(self, term_defs, sets, stable_id):
        for entry in term_defs:
            for term_def in entry.split(";"):
                if term_def:
                    group_name = term_def
                    if not group_name in sets:
                        sets[group_name] = set()
                    if stable_id is not None:
                        sets[group_name].add(stable_id)

    def parse_structure(self):
        parents, children, id_to_descs = self.parse_go_slim()
        (
            parents_protein,
            children_protein,
            id_to_descs_protein,
        ) = self.parse_protein_classes()
        parents.update(parents_protein)
        children.update(children_protein)
        id_to_descs.update(id_to_descs_protein)
        self.structure = (parents, children, id_to_descs)

    def parse_protein_classes(self):
        data_file = os.path.join(self.data_path, "PROTEIN_CLASS_7.0.gz")
        with open_file(data_file, "r") as op:
            parents, children, id_to_descs = _parse_panther_hierarchie(op)
        return parents, children, id_to_descs

    def parse_go_slim(self):
        data_file = os.path.join(self.data_path, "PANTHERGOslim.obo.gz")
        parents = {}
        children = {}
        id_to_descs = {}
        current_id = None
        current_parent = None
        current_name = ""
        with open_file(data_file, "r") as op:
            for line in op:
                line = line.strip()
                if line == "[Term]":
                    if current_id:
                        if not current_parent:
                            current_parent = "top_level"
                        parents[current_id] = current_parent
                        if not current_parent in children:
                            children[current_parent] = set()
                        children[current_parent].add(current_id)
                    id_to_descs[current_id] = current_name
                    current_id = None
                    current_parent = None
                    current_name = ""
                elif line.startswith("id: "):
                    current_id = line[4:].strip()
                elif line.startswith("is_a:"):
                    current_parent = line[len("is_a: ") : line.find(" !")]
                elif line.startswith("name: "):
                    current_name = line[len("name: ") :].strip()
        return parents, children, id_to_descs

    def _set_name_to_group_id(self, set_name):
        protein_classes = re.findall("(PC\d+)", set_name)
        if protein_classes:
            return protein_classes[0]
        go_classes = re.findall("(GO:\d+)", set_name)
        if go_classes:
            return go_classes[0]
        else:
            return "invalid key %s" % set_name

    def get_parent(self, set_name):
        group_id = self._set_name_to_group_id(set_name)
        return self.structure[0][group_id]

    def is_leaf(self, set_name):
        group_id = self._set_name_to_group_id(set_name)
        if not group_id in self.structure[1] or not self.structure[1][group_id]:
            return True
        return False

    def get_url(self, set_name):
        term_id = set_name[: set_name.find("#")]
        return "http://www.pantherdb.org/panther/category.do?categoryAcc=%s" + term_id

    def get_dependencies(self):
        return [
            ppg.FileTimeInvariant(os.path.join(self.data_path, "PTHR7.0_Mouse.bz2")),
            ppg.FileTimeInvariant(os.path.join(self.data_path, "PTHR7.0_Human.bz2")),
            ppg.FunctionInvariant(
                "PantherDBGroups._get_sets", PantherDBDataset._get_sets
            ),
            ppg.FunctionInvariant(
                "PantherDBGroups.parse_terms", PantherDBDataset.parse_terms
            ),
            ppg.FunctionInvariant(
                "PantherDBGroups.parse_structure", PantherDBDataset.parse_structure
            ),
            ppg.FunctionInvariant(
                "_parse_panther_hierarchie", _parse_panther_hierarchie
            ),
            ppg.FunctionInvariant(
                "PantherDBGroups.parse_go_slim", PantherDBDataset.parse_go_slim
            ),
        ]


_david_singelonizer = {}


class _DAVID_Group(_GroupsBase):
    def __new__(cls, name, filename):
        if not name in _david_singelonizer:
            _david_singelonizer[name] = _GroupsBase.__new__(cls)
        return _david_singelonizer[name]

    def __init__(self, name, filename):
        if not hasattr(self, "name"):
            data_dir = os.path.join(
                "/gb/",
                "imt",
                "datasets",
                "DAVID",
                "incoming",
                "6.7",
                "home",
                "davidweb",
                "projects",
                "DAVIDKnowledgebaseRequest",
                "temp",
                "finkernagel_imt.uni-marburg.de",
            )
            full_path = os.path.join(data_dir, filename)
            self.filename = full_path
            _GroupsBase.__init__(self, "DAVID " + name, full_path)
            self.structure = {}

    def get_sets(self, genome):
        if not self.structure:
            self.parse_structure()
        sets = self._parse_term_file()
        res = {}
        for (name, genes) in sets.items():
            if genome.species == "Homo_sapiens":
                genes_human = [x for x in genes if x.startswith("ENSG")]
                genes_mouse = []
            elif genome.species == "Mus_musculus":
                genes_mouse = [x for x in genes if x.startswith("ENSMUSG")]
                genes_human = []
            else:
                # filter to human genes anyw
                genes_human = [x for x in genes if x.startswith("ENSG")]
                genes_mouse = [x for x in genes if x.startswith("ENSMUSG")]

            r = set()
            for g in [
                genes_human,
                genes_mouse,
            ]:  # david has both mouse and human genes in there...
                if g:
                    ok, transformed_genes = self.translate_stable_ids(g, genome)
                    if ok:
                        r.update(transformed_genes)
                    else:
                        raise ValueError(
                            "Could not translate these genes: %s" % (transformed_genes,)
                        )
            res[name] = r
        return res

    @lazy_method
    def _parse_term_file(self):
        res = {}
        op = open_file(self.filename)
        for row in op:
            space_pos = row.find("\t")
            # david_id = int(row[:space_pos])
            ensembl_id = row[:space_pos]
            term = row[space_pos + 1 :].strip()
            if not term in res:
                res[term] = set()
            # res[term].add(david_id)
            res[term].add(ensembl_id)
        return res

    def get_dependencies(self):
        return [
            ppg.FileTimeInvariant(self.filename),
            ppg.FunctionInvariant(
                "_parse_panther_hierarchie", _parse_panther_hierarchie
            ),
        ]

    def get_url(self, name):
        return ""

    def _set_name_to_group_id(self, set_name):
        key = re.findall("((BP|MF)\d+)", set_name)
        if key:
            return key[0][0]
        else:
            return "invalid key"

    def parse_structure(self):
        if self.name == "DAVID PANTHER_BP":
            filename = "/gb/imt/datasets/DAVID/incoming/panther_6.0/Biological_process_6.0.txt.gz"
        elif self.name == "DAVID PANTHER_MF":
            filename = "/gb/imt/datasets/DAVID/incoming/panther_6.0/Molecular_function_6.0.TAB.txt.gz"
        else:
            return
        with open_file(filename, "r") as op:
            self.structure = _parse_panther_hierarchie(op)

    def get_parent(self, set_name):
        if self.structure:
            group_id = self._set_name_to_group_id(set_name)
            return self.structure[0][group_id]

    def is_leaf(self, set_name):
        if self.structure:
            group_id = self._set_name_to_group_id(set_name)
            if not group_id in self.structure[1] or not self.structure[1][group_id]:
                return True
        return False


def DAVID():
    """Get all DAVID defined groups at once"""
    res = []
    for group_name, filename in {
        "GO_TERM_BP_ALL": "ENSEMBL_GENE_ID2GOTERM_BP_ALL.txt",
        "GO_TERM_CC_ALL": "ENSEMBL_GENE_ID2GOTERM_CC_ALL.txt",
        "GO_TERM_MF_ALL": "ENSEMBL_GENE_ID2GOTERM_MF_ALL.txt",
        "KEGG_PATHWAY": "ENSEMBL_GENE_ID2KEGG_PATHWAY.txt",
        "OMIM_DISEASE": "ENSEMBL_GENE_ID2OMIM_DISEASE.txt",
        "PANTHER_BP": "ENSEMBL_GENE_ID2PANTHER_BP_ALL.txt",
        "PANTHER_MF": "ENSEMBL_GENE_ID2PANTHER_MF_ALL.txt",
        "PANTHER_PATHWAY": "ENSEMBL_GENE_ID2PANTHER_PATHWAY.txt",
    }.items():
        res.append(_DAVID_Group(group_name, filename))
    return res


class MacrophageSignatures(_GroupsBase):
    """Xue et al have analyzed 299 microarrays from 29 differently stimulated macrophage conditions and defined 49 coexpression modules - ie
    genes with a high level of correlation, sized 27 to 884 genes. These modules (their eigen-genes...) were correleted to the 29 conditions,
    results in in their suppl. table S2B.
    """

    def __init__(self):
        self.data_path = os.path.join(
            os.path.dirname(__file__), "datasets/XueEtAl2014_Macrophage_Signatures.tsv"
        )
        self.name = "Xue"
        self.raw_sets = None

    def get_dependencies(self):
        return ppg.FileChecksumInvariant(self.data_path)

    def parse(self):
        lookup = {
            "MUTED": "BLOC1S5",
            "KIAA0182": "GSE1",
            "LOC284998": "ENSG00000223806",
            "LOC152195": None,
            "AADACL1": "ENSG00000144959",
            "AFG3L1": "ENSG00000223959",
            "AK3L1": "ENSG00000147853",
            "ALS2CR14": "ENSG00000163596",
            "ALS2CR4": "ENSG00000155755",
            "ANKRD57": "ENSG00000198142",
            "AOF2": "ENSG00000004487",
            "APOB48R": "ENSG00000184730",
            "ARD1A": "ENSG00000268281",
            "ARL17P1": "ENSG00000185829",
            "ARMET": "ENSG00000259880",
            "ARS2": "ENSG00000087087",
            "ATPBD1B": "ENSG00000142751",
            "AXUD1": "ENSG00000144655",
            "AYP1P1": "ENSG00000237659",
            "BAT1": "ENSG00000198563",
            "BAT2D1": "ENSG00000117523",
            "BAT2": "ENSG00000231825",
            "BAT3": "ENSG00000233348",
            "BAT5": "ENSG00000235676",
            "BEXL1": "ENSG00000102409",
            "BHLHB2": "ENSG00000134107",
            "BHLHB3": "ENSG00000123095",
            "BP75": "ENSG00000166164",
            "BRP44": "ENSG00000143158",
            "BRP44L": "ENSG00000060762",
            "BRWD2": "ENSG00000120008",
            "BTBD12": "ENSG00000188827",
            "C10ORF104": "ENSG00000166295",
            "C10ORF119": "ENSG00000197771",
            "C10ORF125": "ENSG00000148803",
            "C10ORF141": "ENSG00000188916",
            "C10ORF26": "ENSG00000166272",
            "C10ORF28": "ENSG00000166024",
            "C10ORF33": "ENSG00000119943",
            "C10ORF46": "ENSG00000151893",
            "C10ORF57": "ENSG00000133678",
            "C10ORF58": "ENSG00000122378",
            "C10ORF59": "ENSG00000184719",
            "C10ORF61": "ENSG00000119977",
            "C10ORF6": "ENSG00000119906",
            "C10ORF78": "ENSG00000156384",
            "C11ORF17": "ENSG00000166452",
            "C11ORF2": "ENSG00000269278",
            "C11ORF46": "ENSG00000152219",
            "C11ORF51": "ENSG00000110200",
            "C11ORF59": "ENSG00000149357",
            "C11ORF60": "ENSG00000269278",
            "C11ORF61": "ENSG00000120458",
            "C11ORF67": "ENSG00000087884",
            "C12ORF11": "ENSG00000064102",
            "C12ORF24": "ENSG00000204856",
            "C12ORF26": "ENSG00000127720",
            "C12ORF30": "ENSG00000111300",
            "C12ORF31": "ENSG00000139233",
            "C12ORF32": "ENSG00000171792",
            "C12ORF35": "ENSG00000174718",
            "C12ORF41": "ENSG00000139620",
            "C12ORF47": "ENSG00000234608",
            "C12ORF59": "ENSG00000165685",
            "C12ORF62": "ENSG00000178449",
            "C13ORF15": "ENSG00000102760",
            "C13ORF18": "ENSG00000102445",
            "C13ORF1": "ENSG00000123178",
            "C13ORF23": "ENSG00000120685",
            "C13ORF27": "ENSG00000151287",
            "C13ORF31": "ENSG00000179630",
            "C13ORF34": "ENSG00000136122",
            "C13ORF7": "ENSG00000152193",
            "C14ORF100": "ENSG00000050130",
            "C14ORF102": "ENSG00000119720",
            "C14ORF104": "ENSG00000165506",
            "C14ORF106": "ENSG00000129534",
            "C14ORF112": "ENSG00000133983",
            "C14ORF118": "ENSG00000089916",
            "C14ORF121": "ENSG00000186648",
            "C14ORF129": "ENSG00000100744",
            "C14ORF131": "ENSG00000022976",
            "C14ORF133": "ENSG00000151445",
            "C14ORF135": "ENSG00000126773",
            "C14ORF138": "ENSG00000100483",
            "C14ORF139": "ENSG00000229645",
            "C14ORF143": "ENSG00000140025",
            "C14ORF147": "ENSG00000165389",
            "C14ORF149": "ENSG00000126790",
            "C14ORF153": "ENSG00000256053",
            "C14ORF156": "ENSG00000119705",
            "C14ORF173": "ENSG00000203485",
            "C14ORF179": "ENSG00000119650",
            "C14ORF21": "ENSG00000196943",
            "C14ORF32": "ENSG00000168175",
            "C14ORF43": "ENSG00000156030",
            "C14ORF4": "ENSG00000119669",
            "C14ORF85": "ENSG00000258730",
            "C15ORF17": "ENSG00000178761",
            "C15ORF21": "ENSG00000179362",
            "C15ORF24": "ENSG00000134153",
            "C15ORF28": None,
            "C15ORF44": "ENSG00000138614",
            "C15ORF63": "ENSG00000242028",
            "C16ORF33": "ENSG00000161981",
            "C16ORF35": "ENSG00000103148",
            "C16ORF42": "ENSG00000007520",
            "C16ORF48": "ENSG00000124074",
            "C16ORF53": "ENSG00000263136",
            "C16ORF56": "ENSG00000102901",
            "C16ORF57": "ENSG00000103005",
            "C16ORF61": "ENSG00000103121",
            "C16ORF63": "ENSG00000133393",
            "C16ORF68": "ENSG00000067365",
            "C16ORF75": "ENSG00000175643",
            "C16ORF79": "ENSG00000182685",
            "C16ORF7": "ENSG00000075399",
            "C17ORF101": "ENSG00000181396",
            "C17ORF39": "ENSG00000141034",
            "C17ORF44": "ENSG00000178977",
            "C17ORF45": "ENSG00000175061",
            "C17ORF48": "ENSG00000170222",
            "C17ORF56": "ENSG00000167302",
            "C17ORF60": "ENSG00000263304",
            "C17ORF63": "ENSG00000173065",
            "C17ORF68": "ENSG00000178971",
            "C17ORF71": "ENSG00000167447",
            "C17ORF81": "ENSG00000170291",
            "C17ORF87": "ENSG00000161929",
            "C17ORF90": "ENSG00000204237",
            "C17ORF91": "ENSG00000186594",
            "C17ORF95": "ENSG00000181038",
            "C18ORF10": "ENSG00000134779",
            "C18ORF19": "ENSG00000177150",
            "C18ORF1": "ENSG00000168675",
            "C18ORF22": "ENSG00000101546",
            "C18ORF55": "ENSG00000075336",
            "C19ORF22": "ENSG00000198858",
            "C19ORF28": "ENSG00000161091",
            "C19ORF29": "ENSG00000226800",
            "C19ORF2": "ENSG00000105176",
            "C19ORF31": None,
            "C19ORF50": "ENSG00000105700",
            "C19ORF56": "ENSG00000105583",
            "C19ORF61": "ENSG00000105771",
            "C19ORF62": "ENSG00000105393",
            "C1ORF107": "ENSG00000117597",
            "C1ORF124": "ENSG00000010072",
            "C1ORF128": "ENSG00000057757",
            "C1ORF144": "ENSG00000055070",
            "C1ORF149": "ENSG00000163875",
            "C1ORF163": "ENSG00000162377",
            "C1ORF166": "ENSG00000090432",
            "C1ORF19": "ENSG00000198860",
            "C1ORF218": "ENSG00000213047",
            "C1ORF24": "ENSG00000135842",
            "C1ORF25": "ENSG00000121486",
            "C1ORF31": "ENSG00000168275",
            "C1ORF38": "ENSG00000130775",
            "C1ORF41": "ENSG00000081870",
            "C1ORF55": "ENSG00000143751",
            "C1ORF57": "ENSG00000135778",
            "C1ORF59": "ENSG00000162639",
            "C1ORF66": "ENSG00000143303",
            "C1ORF69": "ENSG00000181873",
            "C1ORF71": "ENSG00000162852",
            "C1ORF77": "ENSG00000160679",
            "C1ORF83": "ENSG00000116205",
            "C1ORF97": "ENSG00000153363",
            "C20ORF100": "ENSG00000124191",
            "C20ORF107": "ENSG00000213714",
            "C20ORF108": "ENSG00000124098",
            "C20ORF117": "ENSG00000149639",
            "C20ORF11": "ENSG00000101193",
            "C20ORF123": "ENSG00000149635",
            "C20ORF127": "ENSG00000229230",
            "C20ORF160": "ENSG00000101331",
            "C20ORF177": "ENSG00000196227",
            "C20ORF191": "ENSG00000240108",
            "C20ORF20": "ENSG00000101189",
            "C20ORF29": "ENSG00000125843",
            "C20ORF30": "ENSG00000089063",
            "C20ORF3": "ENSG00000101474",
            "C20ORF43": "ENSG00000022277",
            "C20ORF45": "ENSG00000101166",
            "C20ORF4": "ENSG00000131043",
            "C20ORF52": "ENSG00000125995",
            "C20ORF55": "ENSG00000125898",
            "C20ORF7": "ENSG00000101247",
            "C20ORF94": "ENSG00000149346",
            "C21ORF124": "ENSG00000160209",
            "C21ORF24": "ENSG00000223806",
            "C21ORF51": "ENSG00000205670",
            "C21ORF55": "ENSG00000262911",
            "C21ORF57": "ENSG00000182362",
            "C21ORF66": "ENSG00000159086",
            "C21ORF70": "ENSG00000160256",
            "C22ORF13": "ENSG00000138867",
            "C22ORF9": "ENSG00000100364",
            "C2ORF24": "ENSG00000115649",
            "C2ORF25": "ENSG00000168288",
            "C2ORF28": "ENSG00000138085",
            "C2ORF30": "ENSG00000068912",
            "C2ORF32": "ENSG00000119865",
            "C2ORF34": "ENSG00000143919",
            "C2ORF56": "ENSG00000003509",
            "C2ORF64": "ENSG00000183513",
            "C2ORF79": "ENSG00000184924",
            "C2ORF7": "ENSG00000135617",
            "C3ORF10": "ENSG00000254999",
            "C3ORF19": "ENSG00000154781",
            "C3ORF1": "ENSG00000113845",
            "C3ORF21": "ENSG00000173950",
            "C3ORF23": "ENSG00000179152",
            "C3ORF26": "ENSG00000184220",
            "C3ORF31": "ENSG00000264369",
            "C3ORF34": "ENSG00000174007",
            "C3ORF39": "ENSG00000144647",
            "C3ORF59": "ENSG00000180611",
            "C3ORF63": "ENSG00000163946",
            "C3ORF64": "ENSG00000163378",
            "C3ORF75": "ENSG00000163832",
            "C4ORF14": "ENSG00000084092",
            "C4ORF16": "ENSG00000138660",
            "C4ORF18": "ENSG00000164125",
            "C4ORF41": "ENSG00000168538",
            "C5ORF13": "ENSG00000134986",
            "C5ORF21": "ENSG00000113391",
            "C5ORF25": "ENSG00000170085",
            "C5ORF32": "ENSG00000120306",
            "C5ORF35": "ENSG00000155542",
            "C5ORF37": "ENSG00000152359",
            "C5ORF39": "ENSG00000177721",
            "C5ORF41": "ENSG00000164463",
            "C5ORF53": "ENSG00000182700",
            "C5ORF5": "ENSG00000031003",
            "C5ORF62": "ENSG00000256235",
            "C6ORF105": "ENSG00000111863",
            "C6ORF111": "ENSG00000132424",
            "C6ORF115": "ENSG00000146386",
            "C6ORF125": "ENSG00000137288",
            "C6ORF129": "ENSG00000198937",
            "C6ORF130": "ENSG00000124596",
            "C6ORF150": "ENSG00000164430",
            "C6ORF153": "ENSG00000124541",
            "C6ORF160": "ENSG00000203875",
            "C6ORF173": "ENSG00000203760",
            "C6ORF192": "ENSG00000146409",
            "C6ORF59": "ENSG00000203877",
            "C6ORF61": "ENSG00000111877",
            "C6ORF64": "ENSG00000112167",
            "C6ORF66": "ENSG00000123545",
            "C6ORF85": "ENSG00000137266",
            "C6ORF97": "ENSG00000120262",
            "C7ORF11": "ENSG00000168303",
            "C7ORF20": "ENSG00000239857",
            "C7ORF28A": "ENSG00000146574",
            "C7ORF28B": "ENSG00000146574",
            "C7ORF30": "ENSG00000156928",
            "C7ORF36": "ENSG00000241127",
            "C7ORF42": "ENSG00000106609",
            "C7ORF44": "ENSG00000106603",
            "C7ORF47": "ENSG00000160813",
            "C7ORF54": None,
            "C7ORF59": "ENSG00000188186",
            "C7ORF70": "ENSG00000178397",
            "C8ORF38": "ENSG00000156170",
            "C8ORF41": "ENSG00000129696",
            "C8ORF45": "ENSG00000178460",
            "C8ORF55": "ENSG00000263194",
            "C9ORF102": "ENSG00000182150",
            "C9ORF103": "ENSG00000148057",
            "C9ORF109": "ENSG00000231528",
            "C9ORF10OS": "ENSG00000188938",
            "C9ORF119": "ENSG00000175854",
            "C9ORF127": "ENSG00000137103",
            "C9ORF130": "ENSG00000175611",
            "C9ORF164": "ENSG00000187764",
            "C9ORF167": "ENSG00000198113",
            "C9ORF21": "ENSG00000158122",
            "C9ORF23": "ENSG00000164967",
            "C9ORF25": "ENSG00000164970",
            "C9ORF30": "ENSG00000066697",
            "C9ORF46": "ENSG00000107020",
            "C9ORF5": "ENSG00000106771",
            "C9ORF61": "ENSG00000135063",
            "C9ORF6": "ENSG00000119328",
            "C9ORF80": "ENSG00000148153",
            "C9ORF82": "ENSG00000120159",
            "C9ORF86": "ENSG00000196642",
            "C9ORF90": "ENSG00000171169",
            "C9ORF95": "ENSG00000106733",
            "CABC1": "ENSG00000163050",
            "CAMSAP1L1": "ENSG00000118200",
            "CBARA1": "ENSG00000107745",
            "CCDC109A": "ENSG00000156026",
            "CCDC16": "ENSG00000198783",
            "CCDC21": "ENSG00000130695",
            "CCDC45": "ENSG00000263046",
            "CCDC49": "ENSG00000108296",
            "CCDC55": "ENSG00000126653",
            "CCDC56": "ENSG00000183978",
            "CCDC72": "ENSG00000232112",
            "CCDC76": "ENSG00000122435",
            "CCDC99": "ENSG00000040275",
            "CCRK": "ENSG00000156345",
            "CDC2": "ENSG00000138769",
            "CDC2L1": "ENSG00000008128",
            "CDC2L2": "ENSG00000008128",
            "CDC2L5": "ENSG00000065883",
            "CDC2L6": "ENSG00000155111",
            "CDC45L": "ENSG00000093009",
            "CENTA1": "ENSG00000105963",
            "CENTB2": "ENSG00000114331",
            "CENTD2": "ENSG00000186635",
            "CENTG3": "ENSG00000133612",
            "CEP110": "ENSG00000119397",
            "CEP27": "ENSG00000137814",
            "CES4": "ENSG00000228695",
            "CES8": "ENSG00000172824",
            "CHCHD8": "ENSG00000181924",
            "CHCHD9": "ENSG00000186940",
            "CHES1": "ENSG00000108270",
            "CHP": "ENSG00000187446",
            "CICE": "ENSG00000186162",
            "CIP29": "ENSG00000205323",
            "CMAH": "ENSG00000168405",
            "CNO": "ENSG00000186222",
            "COPG": "ENSG00000181789",
            "COX4NB": "ENSG00000131148",
            "CP110": "ENSG00000103540",
            "CRKRS": "ENSG00000167258",
            "CROP": "ENSG00000108848",
            "CRSP2": "ENSG00000180182",
            "CRSP9": "ENSG00000155868",
            "CTGLF3": "ENSG00000204149",
            "CTPS": "ENSG00000171793",
            "CUGBP1": "ENSG00000149187",
            "CUGBP2": "ENSG00000048740",
            "CUTL1": "ENSG00000259938",
            "CXORF12": "ENSG00000268384",
            "CXORF39": "ENSG00000268570",
            "CXORF45": "ENSG00000101901",
            "CYCSL1": "ENSG00000214810",
            "CYLN2": "ENSG00000263226",
            "CYORF15A": "ENSG00000131002",
            "CYORF15B": "ENSG00000131002",
            "CYTSA": "ENSG00000100014",
            "DC36": "ENSG00000196873",
            "DCI": "ENSG00000167969",
            "DDEF2": "ENSG00000151693",
            "DDX39": "ENSG00000123136",
            "DEM1": "ENSG00000164002",
            "DEPDC6": "ENSG00000155792",
            "DJ341D10.1": None,
            "DKFZP434K191": None,
            "DKFZP564O0523": None,
            "DKFZP586I1420": "ENSG00000235859",
            "DKFZP686I15217": "ENSG00000244041",
            "DKFZP686O24166": "ENSG00000188211",
            "DKFZP761E198": "ENSG00000254470",
            "DKFZP761P0423": None,
            "DNCL1": "ENSG00000088986",
            "DULLARD": "ENSG00000175826",
            "DUXAP3": "ENSG00000270552",
            "EBI2": "ENSG00000169508",
            "ECGF1": "ENSG00000025708",
            "ECOP": "ENSG00000154978",
            "EDG1": "ENSG00000170989",
            "EDG4": "ENSG00000064547",
            "EMR4": "ENSG00000268758",
            "FAM108A3": "ENSG00000198658",
            "FAM10A4": "ENSG00000232150",
            "FAM10A7": None,
            "FAM113A": "ENSG00000132635",
            "FAM113B": "ENSG00000179715",
            "FAM116A": "ENSG00000174839",
            "FAM116B": "ENSG00000205593",
            "FAM119A": "ENSG00000144401",
            "FAM119B": "ENSG00000123427",
            "FAM12A": "ENSG00000181562",
            "FAM158A": "ENSG00000100908",
            "FAM164A": "ENSG00000104427",
            "FAM176B": "ENSG00000142694",
            "FAM18B": "ENSG00000171928",
            "FAM23B": "ENSG00000148483",
            "FAM26C": "ENSG00000185933",
            "FAM36A": "ENSG00000203667",
            "FAM38A": "ENSG00000103335",
            "FAM39DP": "ENSG00000146556",
            "FAM39E": "ENSG00000234769",
            "FAM40A": "ENSG00000143093",
            "FAM40B": "ENSG00000128578",
            "FAM44B": "ENSG00000145919",
            "FAM55C": "ENSG00000144815",
            "FAM62B": "ENSG00000117868",
            "FAM86C": "ENSG00000158483",
            "FARSLB": "ENSG00000116120",
            "FBS1": "ENSG00000156860",
            "FBXL10": "ENSG00000272358",
            "FER1L3": "ENSG00000138119",
            "FIS": "ENSG00000203260",
            "FKSG30": "ENSG00000204434",
            "FLJ10081": "ENSG00000114982",
            "FLJ10213": None,
            "FLJ10916": "ENSG00000144115",
            "FLJ10986": "ENSG00000172456",
            "FLJ10996": "ENSG00000125633",
            "FLJ12078": None,
            "FLJ12355": None,
            "FLJ20021": None,
            "FLJ20125": "ENSG00000145723",
            "FLJ20254": None,
            "FLJ20273": None,
            "FLJ20309": None,
            "FLJ20444": "ENSG00000234665",
            "FLJ20489": "ENSG00000211584",
            "FLJ20628": "ENSG00000171103",
            "FLJ20674": None,
            "FLJ20699": "ENSG00000075234",
            "FLJ20718": "ENSG00000155393",
            "FLJ21865": "ENSG00000167280",
            "FLJ21986": "ENSG00000106034",
            "FLJ22639": "ENSG00000225880",
            "FLJ22662": "ENSG00000121316",
            "FLJ23584": "ENSG00000184208",
            "FLJ25006": None,
            "FLJ38482": "ENSG00000170088",
            "FLJ38717": None,
            "FLJ39639": "ENSG00000006194",
            "FLJ39653": "ENSG00000263327",
            "FLJ39827": "ENSG00000184675",
            "FLJ40448": "ENSG00000182376",
            "FLJ40722": "ENSG00000170379",
            "FLJ41423": "ENSG00000255267",
            "FLJ42986": "ENSG00000196460",
            "FLJ44124": None,
            "FLJ45032": "ENSG00000188167",
            "FLJ45244": "ENSG00000235706",
            "FLJ45422": None,
            "FLJ45966": "ENSG00000205959",
            "FLJ46309": None,
            "FLJ46906": "ENSG00000225177",
            "FLJ90757": None,
            "FRAG1": "ENSG00000176208",
            "FRAP1": "ENSG00000198793",
            "FREQ": "ENSG00000107130",
            "FSD1CL": "ENSG00000106701",
            "FTHL11": "ENSG00000237264",
            "FTHL12": "ENSG00000213362",
            "FTHL2": "ENSG00000234975",
            "FTHL3": "ENSG00000223361",
            "FTHL7": "ENSG00000232187",
            "FTHL8": "ENSG00000219507",
            "FVT1": "ENSG00000119537",
            "GALIG": "ENSG00000131981",
            "GARNL4": "ENSG00000132359",
            "GIYD1": "ENSG00000132207",
            "GOLPH4": "ENSG00000173905",
            "GPR109A": "ENSG00000182782",
            "GPR109B": "ENSG00000255398",
            "GPR120": "ENSG00000186188",
            "GPR172A": "ENSG00000185803",
            "GPR175": "ENSG00000163870",
            "GPR177": "ENSG00000116729",
            "GRINL1A": "ENSG00000255529",
            "GSDM1": "ENSG00000167914",
            "GUSBL1": "ENSG00000241549",
            "GVIN1": "ENSG00000254838",
            "HADH2": "ENSG00000268188",
            "HBXIP": "ENSG00000236646",
            "HCCA2": "ENSG00000163374",
            "HCG2P7": None,
            "HDGF2": "ENSG00000166503",
            "HDHD1A": "ENSG00000130021",
            "HIP2": "ENSG00000078140",
            "HISPPD2A": "ENSG00000168781",
            "HLA-A29.1": "ENSG00000224372",
            "HMGB1L1": "ENSG00000124097",
            "HNRPA1L-2": "ENSG00000125820",
            "HNRPA1P4": "ENSG00000206228",
            "HNRPA2B1": "ENSG00000122566",
            "HNRPC": "ENSG00000092199",
            "HNRPH1": "ENSG00000169045",
            "HNRPH3": "ENSG00000096746",
            "HNRPK": "ENSG00000165119",
            "HNRPM": "ENSG00000099783",
            "HNRPR": "ENSG00000125944",
            "HNRPUL1": "ENSG00000105323",
            "HNRPUL2": "ENSG00000214753",
            "HOM-TES-103": "ENSG00000135269",
            "HRASLS3": "ENSG00000176485",
            "HRB": "ENSG00000173744",
            "HSGT1": "ENSG00000122882",
            "HSPC047": None,
            "HSPC111": "ENSG00000048162",
            "HSPC171": "ENSG00000168701",
            "HSPC268": "ENSG00000164898",
            "HSZFP36": None,
            "IFP38": None,
            "IGSF2": "ENSG00000134256",
            "IL1F7": "ENSG00000125571",
            "IL1F9": "ENSG00000136688",
            "IL8RB": "ENSG00000180871",
            "IL8RBP": "ENSG00000229754",
            "IMAA": "ENSG00000260790",
            "ISCA1L": "ENSG00000217416",
            "ISG20L1": "ENSG00000181026",
            "ITGB4BP": "ENSG00000242372",
            "JARID1A": "ENSG00000073614",
            "JARID1D": "ENSG00000012817",
            "JMJD1A": "ENSG00000115548",
            "KIAA0090": "ENSG00000127463",
            "KIAA0174": "ENSG00000182149",
            "KIAA0194": "ENSG00000113716",
            "KIAA0251": "ENSG00000179889",
            "KIAA0261": "ENSG00000062650",
            "KIAA0367": "ENSG00000106772",
            "KIAA0406": "ENSG00000101407",
            "KIAA0415": "ENSG00000242802",
            "KIAA0427": "ENSG00000134030",
            "KIAA0460": "ENSG00000266562",
            "KIAA0495": "ENSG00000227372",
            "KIAA0562": "ENSG00000116198",
            "KIAA0564": "ENSG00000102763",
            "KIAA0776": "ENSG00000014123",
            "KIAA0831": "ENSG00000126775",
            "KIAA0892": "ENSG00000129933",
            "KIAA0913": "ENSG00000214655",
            "KIAA1012": "ENSG00000153339",
            "KIAA1128": "ENSG00000107771",
            "KIAA1160": "ENSG00000240682",
            "KIAA1267": "ENSG00000257382",
            "KIAA1370": "ENSG00000047346",
            "KIAA1530": "ENSG00000163945",
            "KIAA1539": "ENSG00000005238",
            "KIAA1545": "ENSG00000112787",
            "KIAA1600": "ENSG00000151553",
            "KIAA1602": "ENSG00000167566",
            "KIAA1618": "ENSG00000173821",
            "KIAA1632": "ENSG00000152223",
            "KIAA1641": "ENSG00000135976",
            "KIAA1712": "ENSG00000164118",
            "KIAA1751": "ENSG00000142609",
            "KIAA1797": "ENSG00000188352",
            "KIAA1826": "ENSG00000170903",
            "KIAA1862": "ENSG00000133619",
            "KIAA1949": "ENSG00000236428",
            "KIAA2010": "ENSG00000100796",
            "KTELC1": "ENSG00000163389",
            "LASS2": "ENSG00000143418",
            "LASS5": "ENSG00000139624",
            "LASS6": "ENSG00000172292",
            "LAT1-3TM": "ENSG00000260727",
            "LBA1": "ENSG00000168016",
            "LGTN": "ENSG00000265823",
            "LINC01114": "ENSG00000234177",
            "LINCR": "ENSG00000163121",
            "LINS1": "ENSG00000140471",
            "LOC100008589": None,
            "LOC100132288": None,
            "LOC113386": None,
            "LOC123688": "ENSG00000188266",
            "LOC124220": None,
            "LOC124512": None,
            "LOC127295": None,
            "LOC132241": None,
            "LOC133185": None,
            "LOC134997": None,
            "LOC136143": None,
            "LOC143543": None,
            "LOC143666": None,
            "LOC145853": "ENSG00000189227",
            "LOC146517": None,
            "LOC147645": None,
            "LOC147804": None,
            "LOC148915": None,
            "LOC149448": None,
            "LOC151579": None,
            "LOC153561": None,
            "LOC153684": "ENSG00000215068",
            "LOC158160": None,
            "LOC158301": None,
            "LOC158345": None,
            "LOC161527": None,
            "LOC162073": None,
            "LOC196752": None,
            "LOC197135": None,
            "LOC200030": None,
            "LOC201175": None,
            "LOC203547": None,
            "LOC205251": None,
            "LOC220433": None,
            "LOC220686": None,
            "LOC221442": None,
            "LOC221710": None,
            "LOC23117": None,
            "LOC255620": None,
            "LOC255783": None,
            "LOC255809": None,
            "LOC283050": None,
            "LOC283412": None,
            "LOC283874": None,
            "LOC283932": None,
            "LOC284023": "ENSG00000179859",
            "LOC284757": "ENSG00000228340",
            "LOC284821": None,
            "LOC284988": None,
            "LOC285053": None,
            "LOC285074": "ENSG00000231259",
            "LOC285176": None,
            "LOC285550": None,
            "LOC285900": None,
            "LOC286016": None,
            "LOC286208": None,
            "LOC286310": None,
            "LOC338758": "ENSG00000271614",
            "LOC339804": "ENSG00000237651",
            "LOC341457": None,
            "LOC347292": None,
            "LOC347376": None,
            "LOC347544": None,
            "LOC374395": "ENSG00000185475",
            "LOC374443": "ENSG00000256594",
            "LOC387763": None,
            "LOC387820": None,
            "LOC387841": None,
            "LOC387867": None,
            "LOC387934": None,
            "LOC388275": None,
            "LOC388474": None,
            "LOC388524": None,
            "LOC388532": None,
            "LOC388588": None,
            "LOC388621": None,
            "LOC388654": None,
            "LOC388681": None,
            "LOC388720": None,
            "LOC388969": None,
            "LOC389137": None,
            "LOC389203": None,
            "LOC389286": None,
            "LOC389293": None,
            "LOC389435": None,
            "LOC389517": None,
            "LOC389599": None,
            "LOC389641": "ENSG00000246582",
            "LOC389787": None,
            "LOC390354": None,
            "LOC390466": None,
            "LOC391045": None,
            "LOC391656": None,
            "LOC391811": None,
            "LOC399715": "ENSG00000215244",
            "LOC399900": "ENSG00000269089",
            "LOC400027": "ENSG00000273015",
            "LOC400120": None,
            "LOC400455": None,
            "LOC400464": None,
            "LOC400657": "ENSG00000264247",
            "LOC400721": None,
            "LOC400759": None,
            "LOC400890": None,
            "LOC400963": None,
            "LOC400986": None,
            "LOC401019": None,
            "LOC401052": "ENSG00000206567",
            "LOC401115": None,
            "LOC401152": None,
            "LOC401206": None,
            "LOC401218": None,
            "LOC401233": None,
            "LOC401317": "ENSG00000146592",
            "LOC401321": None,
            "LOC401357": "ENSG00000268187",
            "LOC401397": "ENSG00000214194",
            "LOC401622": None,
            "LOC402057": None,
            "LOC402066": None,
            "LOC402221": None,
            "LOC402251": None,
            "LOC402644": None,
            "LOC402694": None,
            "LOC407835": "ENSG00000230626",
            "LOC440055": None,
            "LOC440093": None,
            "LOC440145": "ENSG00000204899",
            "LOC440160": None,
            "LOC440280": None,
            "LOC440341": None,
            "LOC440345": None,
            "LOC440348": "ENSG00000196436",
            "LOC440354": None,
            "LOC440359": None,
            "LOC440589": None,
            "LOC440704": "ENSG00000231175",
            "LOC440731": None,
            "LOC440733": None,
            "LOC440737": None,
            "LOC440900": "ENSG00000234199",
            "LOC440926": None,
            "LOC440927": None,
            "LOC440928": None,
            "LOC440993": "ENSG00000242086",
            "LOC441034": None,
            "LOC441050": None,
            "LOC441087": None,
            "LOC441124": "ENSG00000272645",
            "LOC441126": None,
            "LOC441150": "ENSG00000221821",
            "LOC441155": "ENSG00000058673",
            "LOC441246": None,
            "LOC441377": None,
            "LOC441511": None,
            "LOC441763": None,
            "LOC441775": None,
            "LOC441876": None,
            "LOC442454": None,
            "LOC442535": None,
            "LOC493754": None,
            "LOC54103": "ENSG00000186088",
            "LOC554203": None,
            "LOC606724": None,
            "LOC613037": "ENSG00000185864",
            "LOC641798": None,
            "LOC641825": None,
            "LOC641848": None,
            "LOC641978": None,
            "LOC641992": None,
            "LOC642017": None,
            "LOC642031": "ENSG00000268719",
            "LOC642033": None,
            "LOC642035": None,
            "LOC642082": None,
            "LOC642109": None,
            "LOC642197": None,
            "LOC642210": None,
            "LOC642236": None,
            "LOC642250": None,
            "LOC642255": None,
            "LOC642282": None,
            "LOC642299": None,
            "LOC642333": None,
            "LOC642361": "ENSG00000272447",
            "LOC642468": None,
            "LOC642489": None,
            "LOC642678": None,
            "LOC642732": None,
            "LOC642755": None,
            "LOC642780": None,
            "LOC642817": None,
            "LOC642897": None,
            "LOC642921": None,
            "LOC642934": None,
            "LOC642946": None,
            "LOC642947": None,
            "LOC642989": None,
            "LOC643007": None,
            "LOC643011": None,
            "LOC643031": None,
            "LOC643284": None,
            "LOC643287": None,
            "LOC643300": None,
            "LOC643310": None,
            "LOC643313": None,
            "LOC643357": None,
            "LOC643433": None,
            "LOC643446": None,
            "LOC643452": None,
            "LOC643509": None,
            "LOC643668": None,
            "LOC643790": None,
            "LOC643870": None,
            "LOC643882": None,
            "LOC643913": None,
            "LOC643930": None,
            "LOC643932": None,
            "LOC643949": None,
            "LOC643997": None,
            "LOC644029": None,
            "LOC644033": None,
            "LOC644039": None,
            "LOC644063": None,
            "LOC644090": "ENSG00000242048",
            "LOC644128": None,
            "LOC644131": None,
            "LOC644162": None,
            "LOC644242": None,
            "LOC644250": None,
            "LOC644330": None,
            "LOC644511": None,
            "LOC644590": None,
            "LOC644615": None,
            "LOC644634": "ENSG00000203815",
            "LOC644642": None,
            "LOC644680": None,
            "LOC644739": None,
            "LOC644762": None,
            "LOC644774": None,
            "LOC644799": None,
            "LOC644852": None,
            "LOC644863": None,
            "LOC644869": None,
            "LOC644934": None,
            "LOC644935": None,
            "LOC644979": None,
            "LOC645018": None,
            "LOC645058": None,
            "LOC645138": None,
            "LOC645236": None,
            "LOC645317": None,
            "LOC645362": None,
            "LOC645385": None,
            "LOC645436": None,
            "LOC645466": None,
            "LOC645489": None,
            "LOC645609": None,
            "LOC645683": None,
            "LOC645688": None,
            "LOC645895": None,
            "LOC645899": None,
            "LOC645904": None,
            "LOC645937": None,
            "LOC646064": None,
            "LOC646135": None,
            "LOC646144": None,
            "LOC646195": None,
            "LOC646197": None,
            "LOC646200": None,
            "LOC646463": None,
            "LOC646531": None,
            "LOC646567": None,
            "LOC646630": None,
            "LOC646675": None,
            "LOC646766": None,
            "LOC646786": None,
            "LOC646817": None,
            "LOC646900": None,
            "LOC647000": None,
            "LOC647009": None,
            "LOC647037": None,
            "LOC647041": None,
            "LOC647108": None,
            "LOC647340": None,
            "LOC647346": None,
            "LOC647349": None,
            "LOC647361": None,
            "LOC647389": None,
            "LOC647436": None,
            "LOC647481": None,
            "LOC647673": None,
            "LOC647691": None,
            "LOC647834": None,
            "LOC647856": None,
            "LOC648000": None,
            "LOC648024": None,
            "LOC648099": None,
            "LOC648176": None,
            "LOC648189": None,
            "LOC648210": None,
            "LOC648249": None,
            "LOC648343": None,
            "LOC648366": None,
            "LOC648490": None,
            "LOC648581": None,
            "LOC648605": None,
            "LOC648622": None,
            "LOC648638": None,
            "LOC648695": None,
            "LOC648852": None,
            "LOC648984": None,
            "LOC649044": None,
            "LOC649049": None,
            "LOC649143": None,
            "LOC649150": None,
            "LOC649260": None,
            "LOC649362": None,
            "LOC649422": None,
            "LOC649447": None,
            "LOC649548": None,
            "LOC649555": None,
            "LOC649580": None,
            "LOC649604": None,
            "LOC649661": None,
            "LOC649754": None,
            "LOC649821": None,
            "LOC649841": None,
            "LOC649897": None,
            "LOC649946": None,
            "LOC650029": None,
            "LOC650116": None,
            "LOC650128": None,
            "LOC650152": None,
            "LOC650276": None,
            "LOC650298": None,
            "LOC650369": None,
            "LOC650518": None,
            "LOC650546": None,
            "LOC650557": None,
            "LOC650646": None,
            "LOC650737": None,
            "LOC650832": None,
            "LOC650909": None,
            "LOC650919": None,
            "LOC650950": None,
            "LOC651064": None,
            "LOC651149": None,
            "LOC651202": None,
            "LOC651309": None,
            "LOC651380": None,
            "LOC651436": None,
            "LOC651453": None,
            "LOC651524": None,
            "LOC651575": None,
            "LOC651621": None,
            "LOC651738": None,
            "LOC651816": None,
            "LOC651845": None,
            "LOC651894": None,
            "LOC652071": None,
            "LOC652184": None,
            "LOC652226": None,
            "LOC652264": None,
            "LOC652324": None,
            "LOC652388": None,
            "LOC652455": None,
            "LOC652479": None,
            "LOC652481": None,
            "LOC652489": None,
            "LOC652541": None,
            "LOC652545": None,
            "LOC652595": None,
            "LOC652615": None,
            "LOC652616": None,
            "LOC652624": None,
            "LOC652672": None,
            "LOC652675": None,
            "LOC652685": None,
            "LOC652726": None,
            "LOC652755": None,
            "LOC652844": None,
            "LOC652864": None,
            "LOC652881": None,
            "LOC652968": None,
            "LOC653066": None,
            "LOC653071": None,
            "LOC653073": None,
            "LOC653080": None,
            "LOC653086": None,
            "LOC653103": None,
            "LOC653147": None,
            "LOC653158": None,
            "LOC653171": None,
            "LOC653199": None,
            "LOC653226": None,
            "LOC653232": None,
            "LOC653314": None,
            "LOC653352": None,
            "LOC653354": None,
            "LOC653377": None,
            "LOC653381": None,
            "LOC653382": None,
            "LOC653438": None,
            "LOC653450": None,
            "LOC653468": None,
            "LOC653489": None,
            "LOC653496": None,
            "LOC653505": None,
            "LOC653506": None,
            "LOC653513": None,
            "LOC653566": None,
            "LOC653583": None,
            "LOC653610": None,
            "LOC653631": None,
            "LOC653635": None,
            "LOC653658": None,
            "LOC653702": None,
            "LOC653717": None,
            "LOC653764": None,
            "LOC653773": None,
            "LOC653778": None,
            "LOC653829": None,
            "LOC653874": None,
            "LOC653879": None,
            "LOC653884": None,
            "LOC653888": None,
            "LOC653907": None,
            "LOC653994": None,
            "LOC654000": None,
            "LOC654042": None,
            "LOC654053": None,
            "LOC654069": None,
            "LOC654074": None,
            "LOC654103": None,
            "LOC654121": None,
            "LOC654126": None,
            "LOC654155": None,
            "LOC654174": None,
            "LOC654189": None,
            "LOC654191": None,
            "LOC654194": None,
            "LOC654346": None,
            "LOC727726": None,
            "LOC727759": None,
            "LOC727761": None,
            "LOC727762": None,
            "LOC727820": None,
            "LOC727825": None,
            "LOC727848": None,
            "LOC727935": None,
            "LOC727948": None,
            "LOC728006": None,
            "LOC728037": None,
            "LOC728069": None,
            "LOC728127": None,
            "LOC728153": None,
            "LOC728226": None,
            "LOC728481": None,
            "LOC728492": None,
            "LOC728499": None,
            "LOC728505": None,
            "LOC728518": None,
            "LOC728519": None,
            "LOC728554": "ENSG00000170089",
            "LOC728556": None,
            "LOC728564": None,
            "LOC728565": None,
            "LOC728635": None,
            "LOC728643": None,
            "LOC728653": None,
            "LOC728689": None,
            "LOC728715": "ENSG00000214776",
            "LOC728734": None,
            "LOC728739": None,
            "LOC728744": None,
            "LOC728758": None,
            "LOC728772": None,
            "LOC728888": None,
            "LOC728944": None,
            "LOC728973": None,
            "LOC729008": None,
            "LOC729021": None,
            "LOC729101": None,
            "LOC729148": None,
            "LOC729317": None,
            "LOC729446": None,
            "LOC729466": None,
            "LOC729559": None,
            "LOC729603": "ENSG00000213073",
            "LOC729764": None,
            "LOC729776": None,
            "LOC729843": None,
            "LOC730092": None,
            "LOC730249": None,
            "LOC730316": None,
            "LOC730358": None,
            "LOC730432": None,
            "LOC730455": None,
            "LOC730534": None,
            "LOC730740": None,
            "LOC730744": None,
            "LOC730820": None,
            "LOC730994": None,
            "LOC730995": None,
            "LOC730996": None,
            "LOC731049": None,
            "LOC731096": None,
            "LOC731314": None,
            "LOC731365": None,
            "LOC731486": None,
            "LOC731640": None,
            "LOC731682": None,
            "LOC731777": None,
            "LOC731878": None,
            "LOC731950": None,
            "LOC731985": None,
            "LOC731999": None,
            "LOC732007": None,
            "LOC732075": None,
            "LOC732165": None,
            "LOC732172": None,
            "LOC732371": None,
            "LOC732450": None,
            "LOC81691": "ENSG00000005189",
            "LOC85389": None,
            "LOC85390": None,
            "LOC88523": None,
            "LOC90586": None,
            "LOC90624": None,
            "LOC91561": None,
            "LOC92017": None,
            "LOC93622": "ENSG00000170846",
            "LOC96610": None,
            "LRAP": "ENSG00000164308",
            "LRDD": "ENSG00000177595",
            "LRRC50": "ENSG00000154099",
            "LSM8": "ENSG00000128534",
            "M160": "ENSG00000177675",
            "M6PRBP1": "ENSG00000105355",
            "MAGMAS": "ENSG00000217930",
            "MAK10": "ENSG00000135040",
            "MAP2K1IP1": "ENSG00000109270",
            "MAP3K7IP1": "ENSG00000100324",
            "MAP3K7IP2": "ENSG00000055208",
            "MAPBPIP": "ENSG00000116586",
            "MCART1": "ENSG00000122696",
            "MCMDC1": "ENSG00000111877",
            "METT11D1": "ENSG00000165792",
            "METTL11A": "ENSG00000148335",
            "MGC102966": "ENSG00000214822",
            "MGC10997": None,
            "MGC12760": "ENSG00000215908",
            "MGC15763": "ENSG00000154814",
            "MGC16121": "ENSG00000223749",
            "MGC16169": "ENSG00000145348",
            "MGC16703": None,
            "MGC18216": "ENSG00000140443",
            "MGC27345": "ENSG00000106344",
            "MGC3196": None,
            "MGC33556": "ENSG00000198520",
            "MGC3731": None,
            "MGC40489": None,
            "MGC42367": "ENSG00000196872",
            "MGC4677": "ENSG00000222041",
            "MGC52000": "ENSG00000146556",
            "MGC57346": "ENSG00000204650",
            "MGC70857": "ENSG00000213563",
            "MGC70863": "ENSG00000267966",
            "MGC71993": "ENSG00000219200",
            "MGC72104": None,
            "MGC87042": "ENSG00000105889",
            "MLL4": "ENSG00000105663",
            "MMS19L": "ENSG00000155229",
            "MOBKL2A": "ENSG00000172081",
            "MOBKL2B": "ENSG00000120162",
            "MOBKL2C": "ENSG00000142961",
            "MORG1": "ENSG00000123154",
            "MOSC1": "ENSG00000186205",
            "MRLC2": "ENSG00000118680",
            "MSL3L1": "ENSG00000005302",
            "MTE": "ENSG00000175718",
            "MTMR15": "ENSG00000198690",
            "MTP18": "ENSG00000242114",
            "MUDENG": "ENSG00000053770",
            "MYST1": "ENSG00000103510",
            "MYST3": "ENSG00000083168",
            "NAG18": None,
            "NARG1": "ENSG00000164134",
            "NARG1L": "ENSG00000172766",
            "NAT12": "ENSG00000139977",
            "NAT13": "ENSG00000121579",
            "NAT15": "ENSG00000262621",
            "NAT5": "ENSG00000121579",
            "NAT8B": "ENSG00000144035",
            "NCOA6IP": "ENSG00000137574",
            "NIP30": "ENSG00000172775",
            "NLF2": "ENSG00000205502",
            "NOLA1": "ENSG00000109534",
            "N-PAC": "ENSG00000234799",
            "NPAL3": "ENSG00000001461",
            "NP": "ENSG00000273213",
            "NSBP1": "ENSG00000198157",
            "NSUN5B": "ENSG00000223705",
            "OBFC2A": "ENSG00000173559",
            "OBFC2B": "ENSG00000139579",
            "OKL38": "ENSG00000140961",
            "OR3A4": "ENSG00000180068",
            "ORC2L": "ENSG00000115942",
            "ORC3L": "ENSG00000135336",
            "ORC5L": "ENSG00000164815",
            "ORC6L": "ENSG00000091651",
            "P15RS": "ENSG00000141425",
            "P2RY5": "ENSG00000139679",
            "P704P": None,
            "P76": "ENSG00000125304",
            "P8": "ENSG00000243648",
            "PAPD1": "ENSG00000107951",
            "PCTK3": "ENSG00000117266",
            "PECI": "ENSG00000198721",
            "PFAAP5": "ENSG00000244754",
            "PFTK1": "ENSG00000058091",
            "PGCP": "ENSG00000104324",
            "PHACS": "ENSG00000110455",
            "PHCA": "ENSG00000262970",
            "PIK4CA": "ENSG00000241973",
            "PIP5K2A": "ENSG00000150867",
            "PIP5K2B": "ENSG00000141720",
            "PKM2": "ENSG00000067225",
            "PLDN": "ENSG00000104164",
            "PLEC1": "ENSG00000178209",
            "PLEKHA9": "ENSG00000134297",
            "PMS2L5": "ENSG00000123965",
            "POL3S": "ENSG00000151006",
            "POLS": "ENSG00000106031",
            "PPARBP": "ENSG00000125686",
            "PPIL5": "ENSG00000165501",
            "PPPDE1": "ENSG00000121644",
            "PPPDE2": "ENSG00000100418",
            "PREI3": "ENSG00000115540",
            "PRIC285": "ENSG00000130589",
            "PRKCABP": "ENSG00000100151",
            "PRKCB1": "ENSG00000166501",
            "PRNPIP": "ENSG00000117419",
            "PRO0628": None,
            "PSCD1": "ENSG00000108669",
            "PSCD2": "ENSG00000105443",
            "PSCD4": "ENSG00000100055",
            "PSCDBP": "ENSG00000115165",
            "PXMP3": "ENSG00000164751",
            "RAB7B": None,
            "RABL4": "ENSG00000100360",
            "RAG1AP1": "ENSG00000169241",
            "RAGE": "ENSG00000080823",
            "RAXL1": "ENSG00000173976",
            "RFP": "ENSG00000204713",
            "RG9MTD1": "ENSG00000174173",
            "RICS": "ENSG00000132639",
            "RIPK5": "ENSG00000133059",
            "RNASEN": "ENSG00000113360",
            "RNF160": "ENSG00000198862",
            "RNPC2": "ENSG00000131051",
            "ROD1": "ENSG00000119314",
            "RP11-529I10.4": "ENSG00000115934",
            "RP5-1022P6.2": "ENSG00000215444",
            "RPL13L": "ENSG00000267877",
            "RPRC1": "ENSG00000116871",
            "RPS26L": "ENSG00000196933",
            "RTCD1": "ENSG00000137996",
            "RUNDC2C": "ENSG00000198106",
            "RWDD4A": "ENSG00000182552",
            "SAPS1": "ENSG00000080031",
            "SAPS2": "ENSG00000100239",
            "SAPS3": "ENSG00000110075",
            "SBDSP": "ENSG00000225648",
            "SC4MOL": "ENSG00000052802",
            "SCYL1BP1": "ENSG00000120370",
            "SDCCAG10": "ENSG00000153015",
            "SDCCAG1": "ENSG00000165525",
            "SDHALP1": "ENSG00000185485",
            "SELI": "ENSG00000138018",
            "SELM": "ENSG00000198832",
            "SELO": "ENSG00000073169",
            "SELS": "ENSG00000091490",
            "Sep-02": None,
            "Sep-03": None,
            "Sep-06": None,
            "Sep-07": None,
            "Sep-09": None,
            "Sep-10": "ENSG00000200661",
            "Sep-11": "ENSG00000204571",
            "Sep-15": "ENSG00000261834",
            "SEPT15": None,
            "SEPX1": "ENSG00000198736",
            "SF4": "ENSG00000105705",
            "SFRS10": "ENSG00000136527",
            "SFRS11": "ENSG00000116754",
            "SFRS12": "ENSG00000153914",
            "SFRS13A": "ENSG00000188529",
            "SFRS13B": "ENSG00000154548",
            "SFRS14": "ENSG00000064607",
            "SFRS15": "ENSG00000156304",
            "SFRS16": "ENSG00000104859",
            "SFRS17A": "ENSG00000197976",
            "SFRS1": "ENSG00000164985",
            "SFRS2B": "ENSG00000180771",
            "SFRS2": "ENSG00000161547",
            "SFRS2IP": "ENSG00000139218",
            "SFRS3": "ENSG00000112081",
            "SFRS4": "ENSG00000116350",
            "SFRS5": "ENSG00000100650",
            "SFRS6": "ENSG00000124193",
            "SFRS7": "ENSG00000115875",
            "SFRS8": "ENSG00000061936",
            "SFRS9": "ENSG00000111786",
            "SGK": "ENSG00000118515",
            "SHRM": "ENSG00000138771",
            "SIP1": "ENSG00000065054",
            "SIVA": "ENSG00000184990",
            "SKP1A": "ENSG00000113558",
            "SMA4": "ENSG00000262170",
            "SNHG3-RCC1": "ENSG00000180198",
            "SNX26": "ENSG00000004777",
            "SPG3A": "ENSG00000198513",
            "SPINLW1": "ENSG00000101448",
            "SR140": "ENSG00000163714",
            "SRP14P1": None,
            "SRRM1L": None,
            "STAG3L2": "ENSG00000160828",
            "STS-1": "ENSG00000215946",
            "SUMO1P3": "ENSG00000235082",
            "TADA1L": "ENSG00000152382",
            "TARP": "ENSG00000211689",
            "THEM2": "ENSG00000112304",
            "THOC4": "ENSG00000183684",
            "TIGA1": "ENSG00000224032",
            "TINP1": "ENSG00000164346",
            "TM7SF4": "ENSG00000164935",
            "TMED10P": "ENSG00000254618",
            "TMEM118": "ENSG00000135119",
            "TMEM137": None,
            "TMEM149": "ENSG00000126246",
            "TMEM183B": "ENSG00000163444",
            "TMEM188": "ENSG00000205423",
            "TMEM194": "ENSG00000166881",
            "TMEM1": "ENSG00000160218",
            "TMEM49": "ENSG00000062716",
            "TMEM4": "ENSG00000257727",
            "TMEM77": "ENSG00000156171",
            "TMEM85": "ENSG00000128463",
            "TMEM8": "ENSG00000129925",
            "TMEM93": "ENSG00000127774",
            "TNRC15": "ENSG00000204120",
            "TOP1P2": None,
            "TRA1P2": "ENSG00000203914",
            "TREML3": "ENSG00000184106",
            "TRK1": "ENSG00000171962",
            "TTC15": "ENSG00000171853",
            "TTRAP": "ENSG00000111802",
            "TUBB2C": "ENSG00000188229",
            "TUSC4": "ENSG00000114388",
            "TXNDC14": "ENSG00000213593",
            "TXNDC3": "ENSG00000086288",
            "TXNL2": "ENSG00000108010",
            "U1SNRNPBP": "ENSG00000184209",
            "U2AF1L2": "ENSG00000169249",
            "UBE1C": "ENSG00000144744",
            "UBE1DC1": "ENSG00000081307",
            "UBE1": "ENSG00000130985",
            "UCHL5IP": "ENSG00000268671",
            "UCRC": "ENSG00000184076",
            "UGCGL1": "ENSG00000136731",
            "UGCGL2": "ENSG00000102595",
            "UNC84A": "ENSG00000164828",
            "UNC84B": "ENSG00000100242",
            "UQCRHL": "ENSG00000173660",
            "URG4": "ENSG00000106608",
            "UTX": "ENSG00000147050",
            "VIL2": "ENSG00000092820",
            "VPS24": "ENSG00000249884",
            "VPS26": "ENSG00000122958",
            "WASPIP": "ENSG00000115935",
            "WDR21A": "ENSG00000119599",
            "WDR22": "ENSG00000139990",
            "WDR23": "ENSG00000100897",
            "WDR40A": "ENSG00000198876",
            "WDR51A": "ENSG00000164087",
            "WDR51B": "ENSG00000139323",
            "WDR57": "ENSG00000060688",
            "WDR68": "ENSG00000136485",
            "WDR79": "ENSG00000141499",
            "WDR8": "ENSG00000116213",
            "ZAK": "ENSG00000091436",
            "ZC3H5": "ENSG00000132478",
            "ZNF187": "ENSG00000219891",
            "ZNF285A": "ENSG00000267508",
            "ZNF364": "ENSG00000121848",
            "ZNF509": "ENSG00000168826",
            "ZNF650": "ENSG00000144357",
            "ZNF816A": "ENSG00000180257",
            "ZNF828": "ENSG00000198824",
        }
        filter = set(
            [
                "ENSG00000272645",
                "ENSG00000273213",
                "ENSG00000272358",
                "ENSG00000273015",
                "ENSG00000272447",
                "ENSG00000271614",
                "ENSG00000270552",
            ]
        )
        if not self.raw_sets:
            op = open(self.data_path, "r")
            d = op.readlines()
            in_genes = False
            result = {}
            names_in_order = []
            for ii, line in enumerate(d):
                if line.startswith("Genes"):
                    in_genes = True
                    for name in line.strip().split("\t")[1:]:
                        result[name] = set()
                        names_in_order.append(name)
                elif in_genes:
                    genes_here = line.strip().split("\t")[1:]
                    for gene, name in zip(
                        genes_here, names_in_order[: len(genes_here)]
                    ):
                        gene = gene.strip()
                        if gene and not gene in filter:
                            if gene in lookup:
                                if lookup[gene] and not lookup[gene] in filter:
                                    result[name].add(lookup[gene])
                            else:
                                result[name].add(gene)
            # rename those we have assignments for
            for what, names in [
                ("M1", ["7", "8", "9"]),
                ("M2", ["13", "14", "15"]),
                ("TPP", ["30", "32", "33"]),
            ]:
                for name in names:
                    result[name + " " + what] = result[name]
                    del result[name]
            self.raw_sets = result

    def get_sets(self, genome):
        self.parse()
        result = {}
        any_error = False
        for key in self.raw_sets:
            error, genes = self.translate_hugo_human(
                self.raw_sets[key], genome, ignore_errors=False
            )
            if not error:
                for g in genes:
                    print(g)
                any_error = True
            else:
                result[key] = genes
        if any_error:
            raise ValueError("translation errors")
        return result


class IPA(_GroupsBase):
    def __init__(self):
        self.raw_sets = {}
        self.data_path = os.path.join(os.path.dirname(__file__), "datasets/ipa/")
        self.name = "ipa"

    def get_dependencies(self):
        return [
            # ppg.FileChecksumInvariant(os.path.join(self.data_path, 'Regulator List IPA.txt')),
            ppg.FileChecksumInvariant(
                os.path.join(self.data_path, "Functional Annotation IPA.txt")
            )
        ]

    def parse(self):
        """
        with open(os.path.join(self.data_path, 'Regulator List IPA.txt'), 'r') as op:
            lines = op.read().split("\n")
        lines = [l.split("\t") for l in lines[1:] if l.strip()]
        for l in lines:
            set_name = l[0]
            set_entries = l[2].split(", ")
            self.raw_sets[set_name] = set(set_entries)
        """
        with open(
            os.path.join(self.data_path, "Functional Annotation IPA.txt"), "r"
        ) as op:
            lines = op.read().split("\n")
        lines = [l.split("\t") for l in lines[1:] if l.strip()]
        for l in lines:
            set_name = l[0]
            set_entries = l[1].split(", ")
            self.raw_sets[set_name] = set(set_entries)

    def get_sets(self, genome):
        self.parse()
        result = {}
        any_error = False
        for key in self.raw_sets:
            error, genes_human = self.translate_hugo_human(
                self.raw_sets[key], genome, ignore_errors=True
            )
            if not error:
                for g in genes_human:
                    print(g)
                any_error = True
            else:
                if genome.species != "Homo_sapiens":
                    genes = self.convert_species("Homo_sapiens", genome, genes_human)
                else:
                    genes = genes_human

                result[key] = genes
        if any_error:
            raise ValueError("translation errors")
        return result


class IPA_Regulators(_GroupsBase):
    def __init__(self):
        self.raw_sets = {}
        self.data_path = os.path.join(os.path.dirname(__file__), "datasets/ipa/")
        self.data_file = os.path.join(self.data_path, "Regulator List IPA.txt")
        self.name = "ipa_reg"

    def get_dependencies(self):
        return [
            ppg.FileChecksumInvariant(self.data_file),
            # ppg.FileChecksumInvariant(os.path.join(self.data_path, 'Functional Annotation IPA.txt'))
        ]

    def parse(self):
        """
        with open(os.path.join(self.data_path, 'Regulator List IPA.txt'), 'r') as op:
            lines = op.read().split("\n")
        lines = [l.split("\t") for l in lines[1:] if l.strip()]
        for l in lines:
            set_name = l[0]
            set_entries = l[2].split(", ")
            self.raw_sets[set_name] = set(set_entries)
        """
        with open(self.data_file, "r") as op:
            lines = op.read().split("\n")
        lines = [l.split("\t") for l in lines[1:] if l.strip()]
        for l in lines:
            set_name = l[0]
            set_entries = l[2].split(", ")
            self.raw_sets[set_name] = set(set_entries)

    def get_sets(self, genome):
        self.parse()
        result = {}
        any_error = False
        for key in self.raw_sets:
            error, genes_human = self.translate_hugo_human(
                self.raw_sets[key], genome, ignore_errors=True
            )
            if not error:
                for g in genes_human:
                    print(g)
                any_error = True
            else:
                if genome.species != "Homo_sapiens":
                    genes = self.convert_species("Homo_sapiens", genome, genes_human)
                else:
                    genes = genes_human

                result[key] = genes
        if any_error:
            raise ValueError("translation errors")
        return result


class TFCat(_GroupsBase):
    """interface for the tfcat.ca transcription factor database."""

    def __init__(self, mouse_genome=None):
        self.name = "TFcat"
        self.mouse_genome = mouse_genome
        self.filename = os.path.join(db_dir, "tfcat", "tfcat_export_20170606.tsv")
        self._cache = {}

    def get_dependencies(self):
        return ppg.FileTimeInvariant(self.filename)

    def get_sets(self, genome):
        if genome not in self._cache:
            res = {}
            for name, genes in self.parse(self.filename, genome):
                res[name] = genes
                self._cache[genome] = res
        return self._cache[genome]

    def parse(self, filename, genome):
        import pandas as pd

        if self.mouse_genome is None:
            mouse_genome = ensembl.EnsemblGenome("Mus_musculus", genome.revision)
        else:
            mouse_genome = genome
        print(mouse_genome)
        df = pd.read_csv(filename, sep="\t")
        groups = {}
        for dummy_idx, row in df[["Gene ID", "Judgement"]].iterrows():
            if pd.isnull(row["Gene ID"]):
                continue
            entrez_id = str(int(row["Gene ID"]))  # make sure there's no '.0' at the end
            j = row["Judgement"]
            genes = mouse_genome.entrez_id_to_gene_ids(entrez_id, cache=True)
            if len(genes) == 0:
                genes = mouse_genome.description_to_genes(entrez_id)
                if len(genes) == 0:
                    genes = self.manual_lookup(entrez_id)
            if len(genes):
                if genome.species != "Mus_musculus":
                    genes = self.convert_species("Mus_musculus", genome, genes)
                if j not in groups:
                    groups[j] = set()
                groups[j].update(genes)
        return groups.items()

    @staticmethod
    def manual_lookup(entrez_id):
        manual = {
            11592: None,
            11878: "ENSMUSG00000035277",
            15372: "ENSMUSG00000050100",
            15404: "ENSMUSG00000038236",
            15410: "ENSMUSG00000048763",
            17129: "ENSMUSG00000021540",
            18034: "ENSMUSG00000025225",
            18508: "ENSMUSG00000027168",
            18511: "ENSMUSG00000001497",
            18626: "ENSMUSG00000020893",
            20668: "ENSMUSG00000070643",
            21349: "ENSMUSG00000028717",
            21815: "ENSMUSG00000047407",
            22666: "ENSMUSG00000049672",
            56484: "ENSMUSG00000048756",
            80720: "ENSMUSG00000031860",
            104338: None,
            104385: None,
            104386: None,
            107833: None,
            214162: "ENSMUSG00000002028",
            244349: "ENSMUSG00000031540",
            269424: "ENSMUSG00000025764",
            545848: "ENSMUSG00000076540",
            224829: "ENSMUSG00000064043",
            627232: None,
        }
        if int(entrez_id) in manual:
            res = manual[int(entrez_id)]
            if res is None:
                return set([])
            else:
                return set([res])
        raise KeyError(entrez_id)


FFGroups = GroupsFromDirectory("FF", os.path.join(db_dir, "ff"))
pwc = GMTDataset(
    "PWC", os.path.join(db_dir, "pwc.homo-sapiens-gene-symbol.gmt")
)  # pathway commons


def get_default_groups():
    res = [
        MSigDataset("c2", "v7.2"),
        MSigDataset("c3", "v7.2"),
        MSigDataset("c6", "v7.2"),
        pwc,
        FFGroups,
        IPA(),
    ]
    # if ensembl is not None:
    # res.append(mbf_genomes.EnsemblGO())
    return res
