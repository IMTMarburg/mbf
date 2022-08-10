import pytest
import pandas as pd
import pypipegraph as ppg
from mbf_genomics import genes, DelayedDataFrame
from mbf_genomics.testing import MockGenome
from pypipegraph.testing import force_load
from pathlib import Path


@pytest.mark.usefixtures("new_pipegraph")
class TestDescription:
    def test_simple(self):
        genome = MockGenome(
            pd.DataFrame(
                {
                    "stable_id": ["a", "b", "c"],
                    "chr": "1",
                    "tss": [0, 100, 1000],
                    "tes": [10, 101, 1010],
                }
            ),
            df_genes_meta=pd.DataFrame(
                {
                    "gene_stable_id": ["a", "b", "c"],
                    "description": ["hello", "world", "!"],
                }
            ).set_index("gene_stable_id"),
        )
        g = genes.Genes(genome)
        anno = genes.annotators.Description()
        g += anno
        force_load(g.annotate())
        ppg.run_pipegraph()
        assert "description" in g.df.columns
        assert (
            g.df.sort_values("gene_stable_id")["description"] == ["hello", "world", "!"]
        ).all()

    def test_external_genome(self):
        genome = MockGenome(
            pd.DataFrame(
                {
                    "stable_id": ["a", "b", "c"],
                    "chr": "1",
                    "tss": [0, 100, 1000],
                    "tes": [10, 101, 1010],
                }
            ),
            df_genes_meta=pd.DataFrame(
                {
                    "gene_stable_id": ["a", "b", "c"],
                    "description": ["hello", "world", "!"],
                }
            ).set_index("gene_stable_id"),
        )
        g = DelayedDataFrame("ex", pd.DataFrame({"gene_stable_id": ["a", "c", "b"]}))
        anno = genes.annotators.Description(genome)
        g += anno
        force_load(g.annotate())
        ppg.run_pipegraph()
        assert "description" in g.df.columns
        assert (
            g.df.sort_values("gene_stable_id")["description"] == ["hello", "world", "!"]
        ).all()

    def test_missing_external_genome(self):
        g = DelayedDataFrame("ex", pd.DataFrame({"gene_stable_id": ["a", "c", "b"]}))
        anno = genes.annotators.Description()
        g += anno
        force_load(g.annotate())
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert "ddf had no .genome and no genome was passed to Description" in str(
            g.anno_jobs[anno.get_cache_name()].lfg.exception
        )


@pytest.mark.usefixtures("new_pipegraph")
class TestFromFile:
    def test_simple(self, tmpdir):
        genome = MockGenome(
            pd.DataFrame(
                {
                    "stable_id": ["a", "b", "c"],
                    "chr": "1",
                    "tss": [0, 100, 1000],
                    "tes": [10, 101, 1010],
                }
            ),
            df_genes_meta=pd.DataFrame(
                {
                    "gene_stable_id": ["a", "b", "c"],
                    "description": ["hello", "world", "!"],
                }
            ).set_index("gene_stable_id"),
        )
        g = genes.Genes(genome)
        df_to_add = pd.DataFrame(
            {"testcol": [1, 2, 3], "index_vals": ["a", "b", "d"]}, index=["a", "b", "d"]
        )
        tmp_path = Path(tmpdir) / "dump.tsv"
        df_to_add.to_csv(tmp_path, sep="\t", index=False)
        anno = genes.annotators.FromFile(
            tmp_path,
            columns_to_add=["testcol"],
            index_column_table="index_vals",
            index_column_genes="gene_stable_id",
            fill_value=-1,
        )
        g += anno
        force_load(g.annotate())
        ppg.run_pipegraph()
        print(g.df.index)
        print(g.df)
        assert "testcol" in g.df.columns
        assert g.df.loc[0]["testcol"] == 1
        assert g.df.loc[1]["testcol"] == 2
        assert g.df.loc[2]["testcol"] == -1
        assert len(g.df) == 3
