from mbf.genomics.annotator import Annotator, FromFile
import pandas as pd


class Description(Annotator):
    """Add the description for the genes from genome.

    @genome may be None (default), then the ddf is queried for a '.genome'
    Requires a genome with df_genes_meta - e.g. EnsemblGenomes
    """

    columns = ["description"]

    def __init__(self, genome=None):
        self.genome = genome

    def calc_ddf(self, ddf):
        if self.genome is None:
            try:
                genome = ddf.genome
            except AttributeError:
                raise AttributeError(
                    "ddf had no .genome and no genome was passed to Description"
                )
        else:
            genome = self.genome
        lookup = dict(genome.df_genes_meta["description"].items())
        result = []
        for gene_stable_id in ddf.df["gene_stable_id"]:
            result.append(lookup.get(gene_stable_id, ""))
        return pd.Series(result, index=ddf.df.index)


def GeneStrandedSalmon(*args, **kwargs):
    """Deprecated. use anno_tag_counts.Salmon"""
    raise NotImplementedError("Deprecated. Use anno_tag_counts.Salmon")


# FromFile forwarded to mbf.genomics.annotator.FromFile
FromFile = FromFile


class OldStableIds(Annotator):

    """Add the available old stable ids for each gene
    @genome may be None (default), then the ddf is queried for a '.genome'
    Requires a genome with lookup_stable_id_events. - e.g. EnsemblGenomes

    """

    columns = ["old_gene_stable_ids"]

    def __init__(self, genome=None, max_depth=2):
        self.genome = genome
        self.max_depth = max_depth

    def calc_ddf(self, ddf):
        if self.genome is None:
            try:
                genome = ddf.genome
            except AttributeError:
                raise AttributeError(
                    "ddf had no .genome and no genome was passed to Description"
                )
        else:
            genome = self.genome
        lookup = genome.lookup_stable_id_events.copy()
        # now that's a df {old, new}
        # indexed by old.
        index = {}
        lookup = genome.lookup_stable_id_events
        for old, new in lookup.iterrows():
            for n in new.values:
                for n2 in n:
                    if n2 in index:
                        index[n2].append(old)
                    else:
                        index[n2] = [old]

        def descend(key, accumulator, depth_remaining):
            if depth_remaining == 0:
                return
            try:
                for old in index[key]:
                    accumulator.append(old)
                    descend(old, accumulator, depth_remaining - 1)
            except KeyError:
                return

        result = []
        for gene_stable_id in ddf.df["gene_stable_id"]:
            row = []
            descend(gene_stable_id, row, self.max_depth)
            result.append(", ".join(row))
        return pd.DataFrame({self.columns[0]: result}, index=ddf.df.index)
