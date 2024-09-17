from mbf.genomics.delayeddataframe import DelayedDataFrame
import collections
from ..annotator import Annotator
from ..transcripts import Transcripts
import pandas as pd
import numpy as np


class SummitBase(Annotator):
    pass


class SummitMiddle(SummitBase):
    """Place a summit right in the center (ie. a fake summit"""

    columns = ["summit middle"]
    column_properties = {
        columns[0]: {
            "description": "Fake summit, just the center of the region (given relative to start)"
        }
    }

    def calc(self, df):
        res = []
        for dummy_idx, row in df.iterrows():
            res.append((row["stop"] + row["start"]) / 2 - row["start"])
        return pd.Series(res, dtype=float)


class SummitMACS(SummitBase):
    """Place a summit right from the MACS 'summit' column"""

    def __init__(self):
        self.columns = ["mysummit macs"]
        self.column_properties = {
            self.columns[0]: {
                "description": "Summit position relative to start, as defined by MACS, if not available: fake summit, just the center of the region"
            }
        }
        Annotator.__init__(self)
        pass

    def get_dependencies(self, gr):
        return []

    def calc(self, df):
        # MACS2 stores the summit in absolute coordinates
        abs_summit = df["abs_summit"]
        res = abs_summit - df["start"]
        return pd.DataFrame({self.columns[0]: res})


# from ..genes.anno_tag_counts import GRUnstrandedRust as TagCount
# from ..genes.anno_tag_counts import GRStrandedRust as TagCountStranded
from ..genes.anno_tag_counts import _NormalizationAnno


#class TagCountUnstranded(Annotator):
# that's genes.anno_tag_counts.GRUnstrandedRust you're looking for.

class NormalizationCPM(_NormalizationAnno):
    """Normalize to 1e6 by taking the sum of all regions"""

    def __init__(self, base_column_spec):
        self.name = "CPM(lane)"
        self.normalize_to = 1e6
        super().__init__(base_column_spec)
        self.column_properties = {
            self.columns[0]: {"description": "Tag count normalized to lane tag count"}
        }

    def calc(self, df):
        raw_counts = df[self.raw_column]
        total = max(
            1,
            sum(
                (
                    x.mapped
                    for x in self.raw_anno.aligned_lane.get_bam().get_index_statistics()
                )
            ),
        )
        result = raw_counts * (self.normalize_to / total)
        return pd.Series(result)


class NextTranscript(Annotator):
    """ "Find the closest TSS (any transcript) from the passed genome,
    (or GenomicRegion like (genes/transcripts).filter(....).regions_tss())
    and annotate it's id, gene id, name, and distance"""

    def __init__(self, genome_or_tss_regions):
        if isinstance(genome_or_tss_regions, DelayedDataFrame):
            self.tss = genome_or_tss_regions
        else:
            self.tss = Transcripts(genome_or_tss_regions).regions_tss()
        self.columns = [
            "closest TSS transcript_stable_id",
            "closest TSS gene_stable_id",
            "closest TSS name",
            "closest TSS distance",
            "closest TSS strand",
        ]

    def deps(self, ddf):
        # because deps_annos is not getting the ddf
        return [ddf.anno_jobs[ddf.summit_annotator.get_cache_name()], self.tss.load()]

    def calc_ddf(self, ddf):
        if ddf.genome != self.tss.genome:
            raise ValueError("Mismatched genomes")
        query_positions = ddf.df["start"] + ddf.df[ddf.summit_annotator.columns[0]]
        res = collections.defaultdict(list)
        tss = self.tss
        genome = self.tss.genome
        for chr, qp in zip(ddf.df["chr"], query_positions):
            hit = tss.get_closest_by_start(chr, int(qp))
            if len(hit):
                res["closest TSS transcript_stable_id"].append(
                    hit["transcript_stable_id"].iloc[0]
                )
                res["closest TSS gene_stable_id"].append(hit["gene_stable_id"].iloc[0])
                res["closest TSS name"].append(
                    genome.genes[hit["gene_stable_id"].iloc[0]].name
                )
                res["closest TSS distance"].append(float(hit["start"].iloc[0] - qp))
                res["closest TSS strand"].append(
                    int(
                        hit["tss_direction"]
                        if "tss_direction" in hit
                        else hit["strand"]
                    )
                )
            else:
                res["closest TSS transcript_stable_id"].append("")
                res["closest TSS gene_stable_id"].append("")
                res["closest TSS name"].append("")
                res["closest TSS distance"].append(np.nan)
                res["closest TSS strand"].append(1)

        return pd.DataFrame(res)

    def filter_genes(self, peaks, genes, new_name):
        """Filter a Genes to those that have a peak that has a closest TSS"""
        # implementations considered
        # I can define a filter function, but it's difficult to slot the annotator in place
        # I can define a filter 'definition',  but I obviously don't have the list-of-genes an definition time
        # I can define a filter genes call, taking a name, only drawback is that you might have to filter in multiple steps
        # I can define a genes annotator copying whatever is necessary from, and filter on that ('manually'). (we might need that one anway, someday)
        return genes.filter(
            new_name,
            lambda df: df["gene_stable_id"].isin(peaks.df[self.columns[1]]),
            dependencies=[peaks.add_annotator(self)],
        )
