# flake8: noqa
import pandas as pd
import numpy as np
import pytest
import mbf_genomics.genes as genes
import mbf_genomics.regions as regions
import mbf_genomics.regions.annotators
from mbf_genomics.regions.annotators import SummitMiddle
from .shared import MockGenome, force_load_ddf
import pypipegraph as ppg


@pytest.mark.usefixtures("new_pipegraph")
class TestRegionAnnotationWithGenes:
    def test_anno_next_genes(self):
        genome = MockGenome(
            pd.DataFrame(
                [
                    {
                        "stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5501,
                        "description": "bla",
                        "name": "Fake1",
                    },
                    {
                        "stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                        "name": "Fake2",
                    },
                    {
                        "stable_id": "fake3",
                        "chr": "2",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                        "name": "Fake3",
                    },
                ]
            ),
            df_transcripts=pd.DataFrame(
                [
                    {
                        "transcript_stable_id": "fake1a",
                        "gene_stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5000,
                        "tes": 5500,
                        "description": "bla",
                        "exons": [(5000, 5500)],
                        "name": "Fake1",
                    },
                    {
                        "transcript_stable_id": "fake1b",
                        "gene_stable_id": "fake1",
                        "chr": "1",
                        "strand": 1,
                        "tss": 5500,
                        "tes": 5501,
                        "description": "bla",
                        "exons": [(5500, 5501)],
                        "name": "Fake1b",
                    },
                    {
                        "transcript_stable_id": "fake2a",
                        "gene_stable_id": "fake2",
                        "chr": "1",
                        "strand": -1,
                        "tss": 5400,
                        "tes": 4900,
                        "description": "bla",
                        "exons": [(4900, 5400)],
                        "name": "Fake2a",
                    },
                ]
            ),
        )

        def sample_data():
            df = pd.DataFrame(
                {
                    "chr": ["1", "2", "1", "1", "3", "5"],
                    "start": [10, 100, 6000, 5400, 10000, 100000],
                    "stop": [12, 110, 6110, 5500, 11110, 111110],
                }
            )
            # df = df.assign(summit=(df["stop"] - df["start"]) // 2)
            return df

        a = regions.GenomicRegions("shu", sample_data, [], genome)
        # a.summit_annotator = SummitMiddle()
        anno = mbf_genomics.regions.annotators.NextTranscript(genome)
        a.add_annotator(anno)
        force_load_ddf(a)

        my_genes = genes.Genes(genome)
        fg = anno.filter_genes(a, my_genes, 'filtered_genes')
        force_load_ddf(fg)
        ppg.run_pipegraph()

        # remember, these are sorted chr, start
        assert (
            a.df["closest TSS name"]
            == np.array(["Fake1", "Fake2", "Fake1", "", "", ""])
        ).all()
        should = [
            1.0 * (5000 - (11)),
            -1.0 * (5450 - 5400),
            -1.0 * ((6000 + 6110) / 2 - 5500),
            np.nan,
            np.nan,
            np.nan,
        ]
        assert (
            (a.df["closest TSS distance"] == should)
            | np.isnan(a.df["closest TSS distance"])
        ).all()

        print(dir(fg))
        assert len(fg.df) == 2
        assert (fg.df['gene_stable_id'] == [ 'fake2','fake1']).all() # fake1 gene is earlier with it's TES than fake1
