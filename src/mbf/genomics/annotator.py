from abc import ABC
from typing import List
import hashlib
import pandas as pd
import pypipegraph as ppg
from .util import freeze
from pathlib import Path
import numpy as np

annotator_singletons = {"lookup": []}


class Annotator(ABC):
    def __new__(cls, *args, **kwargs):
        cn = cls.__name__
        if ppg.inside_ppg():
            if not hasattr(ppg.util.global_pipegraph, "_annotator_singleton_dict"):
                ppg.util.global_pipegraph._annotator_singleton_dict = {"lookup": []}
            singleton_dict = ppg.util.global_pipegraph._annotator_singleton_dict
        else:
            singleton_dict = annotator_singletons
        if not cn in singleton_dict:
            singleton_dict[cn] = {}
        key = cls._prep_init_params_for_freezing(*args, **kwargs)
        for k, v in key.items():
            try:
                key[k] = freeze(v)
            except Exception as e:
                raise ValueError(e, "key was", k, "value is", v)
        key = tuple(sorted(key.items()))
        if not key in singleton_dict[cn]:
            singleton_dict[cn][key] = object.__new__(cls)
            singleton_dict["lookup"].append(singleton_dict[cn][key])

        res = singleton_dict[cn][key]
        if res.__class__ != cls:
            raise ValueError(
                "Two Annotator classes with same class name will lead to adding the wrong annotator. Rename one of them"
            )
        return res

    @classmethod
    def _prep_init_params_for_freezing(cls, *args, **kwargs):
        """collect what we need to freeze in order to singletonize this Annotator"""
        key = {}
        for ii in range(0, len(args)):
            key["arg_%i" % ii] = args[ii]
        key.update(kwargs)
        return key

    def __hash__(self):
        return hash(self.get_cache_name())

    def __str__(self):
        return "Annotator %s" % self.columns[0]

    def __repr__(self):
        return "Annotator(%s)" % self.columns[0]

    def __freeze__(self):
        return "Annotator(%s)" % self.columns[0]

    def get_cache_name(self):
        if hasattr(self, "cache_name"):
            return self.cache_name
        else:
            return self.columns[0]

    def calc(self, df):
        raise NotImplementedError()  # pragma: no cover

    def deps(self, ddf):
        """Return ppg.jobs"""
        return []

    def dep_annos(self):
        """Return other annotators"""
        return []


class Constant(Annotator):
    def __init__(self, column_name, value):
        self.columns = [column_name]
        self.value = value

    def calc(self, df):
        return pd.DataFrame({self.columns[0]: self.value}, index=df.index)


class FromFile(Annotator):
    def __init__(
        self,
        tablepath: Path,
        columns_to_add: List[str],
        index_column_table: str = "gene_stable_id",
        index_column_genes: str = "gene_stable_id",
        fill_value: float = None,
        is_tsv: bool = False,
    ):
        """
        Adds arbitrary columns from a table.

        This requires that both the table and the ddf have a common column on
        which we can index.

        Parameters
        ----------
        tablepath : Path
            Path to table with additional columns.
        columns_to_add : List[str]
            List of columns to append.
        index_column_table : str, optional
            Index column in table, by default "gene_stable_id".
        index_column_genes : str, optional
            Index column in ddf to append to, by default "gene_stable_id".
        fill_value : float, optonal
            Value to fill for missing rows, defaults to np.NaN.
        is_tsv : bool
            If the input file is a .tsv file regardless of the suffix.
        """
        self.tablepath = tablepath
        self.columns = columns_to_add
        self.index_column_table = index_column_table
        self.index_column_genes = index_column_genes
        self.fill = fill_value if fill_value is not None else np.NaN
        self.is_tsv = is_tsv

    def parse(self):
        if (
            (self.tablepath.suffix == ".xls") or (self.tablepath.suffix == ".xlsx")
        ) and not self.is_tsv:
            return pd.read_excel(self.tablepath)
        else:
            return pd.read_csv(self.tablepath, sep="\t")

    def get_cache_name(self):
        suffix = f"{self.tablepath.name}_{self.columns[0]}".encode("utf-8")
        return f"FromFile_{hashlib.md5(suffix).hexdigest()}"

    def calc_ddf(self, ddf):
        """Calculates the ddf to append."""
        df_copy = ddf.df.copy()
        if self.index_column_genes not in df_copy.columns:
            raise ValueError(
                f"Column {self.index_column_genes} not found in ddf index, found was:\n{[str(x) for x in df_copy.columns]}."
            )
        df_in = self.parse()
        if self.index_column_table not in df_in.columns:
            raise ValueError(
                f"Column {self.index_column_table} not found in table, found was:\n{[str(x) for x in df_in.columns]}."
            )
        for column in self.columns:
            if column not in df_in.columns:
                raise ValueError(
                    f"Column {column} not found in table, found was:\n{[str(x) for x in df_in.columns]}."
                )
        df_copy.index = df_copy[self.index_column_genes]
        df_in.index = df_in[self.index_column_table]
        df_in = df_in.reindex(df_copy.index, fill_value=self.fill)
        df_in = df_in[self.columns]
        df_in.index = ddf.df.index
        return df_in

    def deps(self, ddf):
        """Return ppg.jobs"""
        return ppg.FileInvariant(self.tablepath)
