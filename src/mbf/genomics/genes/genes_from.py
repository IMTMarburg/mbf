import numpy as np
import pypipegraph as ppg
from .genes import Genes
from ..util import read_pandas, freeze


def FromDifference(name, a, b, sheet_name="Differences"):

    """a minus b"""

    def do_load(df):
        remove_ids = set(b.df["gene_stable_id"])
        keep = ~np.array(
            [stable_id in remove_ids for stable_id in a.df["gene_stable_id"]],
            dtype=bool,
        )
        return keep

    if a.load_strategy.build_deps:
        deps = [
            a.load(),
            b.load(),
            ppg.ParameterInvariant(
                "Genes_%s_parents" % name, (a.name, b.name)
            ),  # so if you swap out the gr, it's detected...
        ]
    else:
        deps = []

    res = a.filter(
        name,
        do_load,
        dependencies=deps,
        sheet_name=sheet_name,
        vid=["Difference", a.vid, b.vid],
    )
    res.parent = a
    return res


def FromIntersection(name, gene_sets, sheet_name="Intersections"):
    if not isinstance(sheet_name, str):
        raise ValueError(
            "sheet name must be string - note that FromIntersection"
            "takes a list of genes as 2nd parameter"
        )
    parents = []
    vid = _vid_from_genes_sets(gene_sets, "Intersection")
    for g in gene_sets:
        parent = g.parent
        while parent.parent:
            parent = parent.parent
        parents.append(parent)
    parent_names = set([x.name for x in parents])
    if len(parent_names) != 1:  # pragma: no cover
        raise ValueError(
            "Trying to combine gene lists from different parents - currently not supported"
        )

    def in_set(df):
        accepted_ids = set.intersection(
            *[set(g.df["gene_stable_id"]) for g in gene_sets]
        )
        return np.array([x in accepted_ids for x in df["gene_stable_id"]], dtype=bool)

    return parents[0].filter(
        name,
        in_set,
        dependencies=[g.load() for g in gene_sets],
        sheet_name=sheet_name,
        vid=vid,
    )


def _vid_from_genes_sets(genes_sets, prefix="filtered", sheet_name=None):
    if sheet_name is not None:
        vid = [prefix, sheet_name]
    else:
        vid = [prefix]
    for g in genes_sets:
        if g.vid is not None:
            if isinstance(g.vid, str):
                vid.append(g.vid)
            else:
                vid = vid + g.vid
    return vid


def _from_filtered_genes(name, do_filter, genes_sets, sheet_name, vid):
    if genes_sets[0].load_strategy.build_deps:
        deps = ([o.load() for o in genes_sets],)
    else:
        deps = []
    return Genes(genes_sets[0].genome).filter(
        name, do_filter, dependencies=deps, sheet_name=sheet_name, vid=vid
    )


def FromAny(name, genes_sets, sheet_name=None):
    def do_filter(df):
        seen = set()
        for o in genes_sets:
            seen.update(o.df["gene_stable_id"])
        return np.array([stable_id in seen for stable_id in df["gene_stable_id"]])

    return _from_filtered_genes(
        name,
        do_filter,
        genes_sets,
        sheet_name,
        vid=_vid_from_genes_sets(genes_sets, "filtered", sheet_name),
    )


def FromAll(name, genes_sets, sheet_name=None):
    # def filter_to_those_occuring_in_all_filtered_sets(
    def do_filter(df):
        ok = set.intersection(
            *[set(o.df["gene_stable_id"].unique()) for o in genes_sets]
        )
        return np.array([stable_id in ok for stable_id in df["gene_stable_id"]])

    return _from_filtered_genes(
        name,
        do_filter,
        genes_sets,
        sheet_name,
        vid=_vid_from_genes_sets(genes_sets, "filtered", sheet_name),
    )


def FromNone(name, genes_sets, sheet_name=None):
    # def filter_to_those_not_occuring_in_any_filtered_sets(
    def do_filter(df):
        seen = set()
        for o in genes_sets:
            seen.update(o.df["gene_stable_id"])
        return np.array(
            ~np.array(
                [stable_id in seen for stable_id in df["gene_stable_id"]], dtype=bool
            )
        )

    return _from_filtered_genes(
        name,
        do_filter,
        genes_sets,
        sheet_name,
        vid=_vid_from_genes_sets(genes_sets, "filtered", sheet_name),
    )


def FromFile(
    name, genome, table_filename, column_name="gene_stable_id", sheet_name=None
):
    """Filter Genes(genome) to those occuring in the table_filename"""

    table_filename = str(table_filename)

    def filter(genes_df):

        df = read_pandas(table_filename)
        seen = df[column_name].values
        return np.array(
            [str(x) in seen for x in genes_df["gene_stable_id"]], dtype=bool
        )

    g = Genes(genome)
    if g.load_strategy.build_deps:
        deps = [ppg.FileChecksumInvariant(table_filename)]
    else:
        deps = []

    return g.filter(name, filter, dependencies=deps, sheet_name=sheet_name)


def FromFileOfTranscripts(
    name, genome, table_filename, column_name="transcript_stable_id"
):
    """Filter Genes(genome) to those whose transcripts occur in the table_filename"""

    def filter(genes_df):
        df = read_pandas(table_filename)
        transcripts = df[column_name]
        seen = set()
        for transcript_stable_id in transcripts:
            seen.add(genome.transcripts[transcript_stable_id].gene_stable_id)

        return np.array([x in seen for x in genes_df["gene_stable_id"]], dtype=bool)

    g = Genes(genome)
    if g.load_strategy.build_deps:
        deps = [ppg.FileChecksumInvariant(table_filename)]
    else:
        deps = []
    return g.filter(name, filter, dependencies=deps)


def FromBiotypes(genome, allowed_biotypes):
    def filter(genes_df):
        ok = np.zeros(genes_df.shape[0], dtype=bool)
        for x in allowed_biotypes:
            ok = ok | (genes_df["biotype"] == x)
        return ok

    return Genes(genome).filter(
        "%s_with_%s" % (genome.name, ",".join(allowed_biotypes)), filter
    )


def FromNames(
    name,
    genome,
    list_or_callback,
    sheet_name=None,
    vid=None,
    manual_lookup=None,
    ignore_unmatched=False,
):
    """Filter Genes(genome) to those occuring in the list (or the list
    returned from the callback.

    If manual_lookup is set, genes will be first replaced with those before lookup
    in the genome. Entries resolving to none will be excluded.
    Manual lookup may contain unused entries

    """

    def filter(genes_df, manual_lookup=manual_lookup):
        if hasattr(list_or_callback, "__call__"):
            ll = list_or_callback()
        else:
            ll = list_or_callback
        if manual_lookup is not None:
            manual_lookup = {
                k.upper(): v.upper() if v is not None else None
                for (k, v) in manual_lookup.items()
            }
            ll = [manual_lookup.get(x.upper(), x.upper()) for x in ll]
            ll = [x for x in ll if x is not None]

        stable_ids = set(genes_df["gene_stable_id"])
        ids_seen = [x for x in ll if x in stable_ids]
        names = [x for x in ll if x not in stable_ids]
        seen = set([x.upper() for x in names])
        unused = seen.difference(genes_df["name"].str.upper())
        if unused and not ignore_unmatched:
            raise ValueError(
                "the following gene names were not found in the genome:\n"
                + "\n".join(sorted(unused))
            )
        res = np.zeros(len(genes_df), bool)
        res |= genes_df["gene_stable_id"].isin(ids_seen)
        res |= genes_df["name"].str.upper().isin(seen)
        return res

    g = Genes(genome)
    if g.load_strategy.build_deps:
        if hasattr(list_or_callback, "__call__"):
            deps = [ppg.FunctionInvariant(name + "_list_or_callback", list_or_callback)]
        else:
            deps = [
                ppg.ParameterInvariant(name + "_list_or_callback", list_or_callback)
            ]
        deps.append(
            ppg.ParameterInvariant(name + "_manual_lookup", freeze(manual_lookup))
        )
    else:
        deps = []

    return g.filter(name, filter, dependencies=deps, sheet_name=sheet_name, vid=vid)
