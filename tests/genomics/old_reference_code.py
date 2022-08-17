from mbf.genomics.annotator import Annotator
import hashlib
import numpy as np
import pandas as pd


class NormalizationCPMBiotypes(Annotator):
    """Tormalize to 1e6 by taking the sum of all [biotype, biotype2] genes.
    All other genes receive nan as their normalized value"""

    def __init__(self, raw_anno, biotypes):
        self.genome = raw_anno.genome
        if not isinstance(biotypes, tuple):  # pragma: no cover
            raise ValueError("biotypes must be a tuple")
        self.biotypes = biotypes
        self.raw_anno = raw_anno
        self.vid = raw_anno.vid
        self.normalize_to = 1e6
        self.aligned_lane = raw_anno.aligned_lane
        self.columns = [
            self.raw_anno.columns[0] + " CPM(%s)" % (", ".join(sorted(biotypes)))
        ]
        self.cache_name = hashlib.md5(self.columns[0].encode("utf-8")).hexdigest()
        self.column_properties = {
            self.columns[0]: {
                "description": "Tag count inside protein coding (all if no protein coding transcripts) exons, normalized to 1e6 across genes in biotypes %s"
                % (biotypes,)
            }
        }

    def dep_annos(self):
        return [self.raw_anno]

    def calc(self, df):
        raw_counts = df[self.raw_anno.columns[0]].copy()
        ok = np.zeros(len(df), bool)
        for biotype in self.biotypes:
            ok |= df["biotype"] == biotype
        raw_counts[~ok] = np.nan
        total = float(raw_counts[ok].sum())
        result = raw_counts * (self.normalize_to / total)
        return pd.DataFrame({self.columns[0]: result})


class NormalizationTPMBiotypes(Annotator):
    """TPM, but only consider genes matching one of the biotypes
    All other genes receive nan as their normalized value"""

    def __init__(self, raw_anno, biotypes):
        self.genome = raw_anno.genome
        if not isinstance(biotypes, tuple):  # pragma: no cover
            raise ValueError("biotypes must be a tuple")
        self.biotypes = biotypes
        self.raw_anno = raw_anno
        self.vid = raw_anno.vid
        self.normalize_to = 1e6
        self.aligned_lane = raw_anno.aligned_lane
        self.columns = [
            self.raw_anno.columns[0] + " tpm(%s)" % (", ".join(sorted(biotypes)))
        ]
        self.cache_name = hashlib.md5(self.columns[0].encode("utf-8")).hexdigest()
        self.column_properties = {
            self.columns[0]: {
                "description": "transcripts per million, considering only biotypes %s"
                % (biotypes,)
            }
        }

    def dep_annos(self):
        return [self.raw_anno]

    def calc(self, df):
        raw_counts = df[self.raw_anno.columns[0]].copy()
        ok = np.zeros(len(df), bool)
        for biotype in self.biotypes:
            ok |= df["biotype"] == biotype
        raw_counts[~ok] = np.nan

        length_by_gene = self.raw_anno.interval_strategy.get_interval_lengths_by_gene(
            self.genome
        )
        result = np.zeros(raw_counts.shape, float)
        for ii, gene_stable_id in enumerate(df["gene_stable_id"]):
            result[ii] = raw_counts[ii] / float(length_by_gene[gene_stable_id])
        total = float(result[ok].sum())  # result.sum would be nan!
        factor = self.normalize_to / total
        result = result * factor
        return pd.DataFrame({self.columns[0]: result})
