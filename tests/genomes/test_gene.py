import pytest
import numpy as np
from mbf_genomes import HardCodedGenome
import pandas as pd

default_chr_lengths = {
    "1": 100_000,
    "2": 200_000,
    "3": 300_000,
    "4": 400_000,
    "5": 500_000,
}


def DummyGenome(df_genes, df_transcripts=None):

    df_genes = df_genes.rename(columns={"stable_id": "gene_stable_id"})
    if not "start" in df_genes.columns:
        starts = []
        stops = []
        for idx, row in df_genes.iterrows():
            if row["strand"] == 1:
                starts.append(row["tss"])
                stops.append(row["tes"])
            else:
                starts.append(row["tes"])
                stops.append(row["tss"])
        df_genes = df_genes.assign(start=starts, stop=stops)
    if not "biotype" in df_genes.columns:
        df_genes = df_genes.assign(biotype="protein_coding")
    if not "name" in df_genes.columns:
        df_genes = df_genes.assign(name=df_genes.index)
    df_genes = df_genes.sort_values(["chr", "start"])
    df_genes = df_genes.set_index("gene_stable_id")
    if df_transcripts is not None:
        if not "biotype" in df_transcripts.columns:
            df_transcripts = df_transcripts.assign(biotype="protein_coding")
        if not "name" in df_transcripts.columns:
            df_transcripts = df_transcripts.assign(name=df_transcripts.index)
        if "exons" in df_transcripts.columns:
            if len(df_transcripts["exons"].iloc[0]) == 3:
                df_transcripts = df_transcripts.assign(
                    exons=[(x[0], x[1]) for x in df_transcripts["exons"]]
                )
            df_transcripts = df_transcripts.assign(
                exon_stable_ids=[
                    "exon_%s_%i" % (idx, ii)
                    for (ii, idx) in enumerate(df_transcripts["exons"])
                ]
            )
        df_transcripts = df_transcripts.set_index("transcript_stable_id")
    return HardCodedGenome("dummy", default_chr_lengths, df_genes, df_transcripts, None)


def test_transcript_get_introns():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4900, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3750, 4000)],
                    [(4900, 5000), (5100, 5400)],
                    [(4900, 5000), (5100, 5200), (5222, 5400)],
                ],
            }
        ),
    )
    assert genome.transcripts["trans1a"].introns == [(3000, 3100)]
    assert genome.transcripts["trans1b"].introns == [(3500, 3750), (4000, 4900)]
    assert genome.transcripts["trans2"].introns == [(5000, 5100)]
    assert genome.transcripts["trans3"].introns == [(5000, 5100), (5200, 5222)]


def test_gene_get_introns_merging():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4900, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4800)],
                    [(3000, 3500), (3750, 4000)],
                    [(4900, 5000), (5100, 5400)],
                    [(4900, 5000), (5100, 5200), (5222, 5400)],
                ],
            }
        ),
    )
    print(genome.genes["fake1"].introns_strict)

    def as_list(introns):
        return [list(introns[0]), list(introns[1])]

    assert as_list(genome.genes["fake1"].introns_strict) == [[4800], [4900]]
    # assert genome.transcripts["trans1b"].introns == [(3500, 3750), (4000, 4900)]
    assert as_list(genome.genes["fake2"].introns_strict) == [[5000], [5100]]
    assert as_list(genome.genes["fake3"].introns_strict) == [[5000, 5200], [5100, 5222]]
    assert as_list(genome.genes["fake1"].introns_all) == [
        [3000, 4800, 3500, 4000],
        [3100, 4900, 3750, 4900],
    ]


def test_intronify_more_complex():
    transcript = {
        "chr": "2R",
        "exons": [
            (14_243_005, 14_244_766),
            (14_177_040, 14_177_355),
            (14_176_065, 14_176_232),
            (14_175_632, 14_175_893),
            (14_172_742, 14_175_244),
            (14_172_109, 14_172_226),
            (14_170_836, 14_172_015),
            (14_169_750, 14_170_749),
            (14_169_470, 14_169_683),
            (14_169_134, 14_169_402),
            (14_167_751, 14_169_018),
            (14_166_570, 14_167_681),
        ],
        "gene_stable_id": "FBgn0010575",
        "start": 14_166_570,
        "stop": 14_244_766,
        "strand": -1,
        "transcript_stable_id": "FBtr0301547",
    }
    gene = {
        "biotype": "protein_coding",
        "chr": "2R",
        "description": "CG5580 [Source:FlyBase;GeneId:FBgn0010575]",
        "name": "sbb",
        "stable_id": "FBgn0010575",
        "strand": -1,
        "tes": 14_166_570,
        "tss": 14_244_766,
    }
    genome = DummyGenome(pd.DataFrame([gene]), pd.DataFrame([transcript]))
    g = genome.transcripts["FBtr0301547"]
    assert g.gene_stable_id == "FBgn0010575"
    introns = g.introns
    assert (
        np.array(introns)
        == [
            (14_167_681, 14_167_751),
            (14_169_018, 14_169_134),
            (14_169_402, 14_169_470),
            (14_169_683, 14_169_750),
            (14_170_749, 14_170_836),
            (14_172_015, 14_172_109),
            (14_172_226, 14_172_742),
            (14_175_244, 14_175_632),
            (14_175_893, 14_176_065),
            (14_176_232, 14_177_040),
            (14_177_355, 14_243_005),
        ]
    ).all()
    assert genome.get_chromosome_lengths() == default_chr_lengths
    assert genome.job_genes() is None


def test_intron_intervals_raises_on_inverted():
    transcript = {
        "chr": "2R",
        "exons": [
            (14_243_005, 14_244_766),
            (14_177_040, 14_177_355),
            (14_176_065, 14_176_232),
            (14_175_632, 14_175_893),
            (14_172_742, 14_175_244),
            (14_172_109, 14_172_226),
            (14_172_015, 14_170_836),  # inverted
            (14_169_750, 14_170_749),
            (14_169_470, 14_169_683),
            (14_169_134, 14_169_402),
            (14_167_751, 14_169_018),
            (14_166_570, 14_167_681),
        ],
        "gene_stable_id": "FBgn0010575",
        "start": 14_166_570,
        "stop": 14_244_766,
        "strand": -1,
        "transcript_stable_id": "FBtr0301547",
    }
    gene = {
        "biotype": "protein_coding",
        "chr": "2R",
        "description": "CG5580 [Source:FlyBase;GeneId:FBgn0010575]",
        "name": "sbb",
        "stable_id": "FBgn0010575",
        "strand": -1,
        "tes": 14_166_570,
        "tss": 14_244_766,
    }
    genome = DummyGenome(pd.DataFrame([gene]), pd.DataFrame([transcript]))
    g = genome.transcripts["FBtr0301547"]
    with pytest.raises(ValueError):
        g.introns


def test_get_gene_introns():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5500,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4910, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    one = genome.genes["fake1"].introns_strict
    assert len(one[0]) == 0

    two = genome.genes["fake2"].introns_strict

    def as_list(introns):
        return [list(introns[0]), list(introns[1])]

    assert as_list(two) == [[4900, 5000, 5400], [4910, 5100, 5500]]


def test_get_gene_exons_merged():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
                {
                    "stable_id": "fake4",
                    "chr": "2",
                    "strand": -1,
                    "tss": 6400,
                    "tes": 5900,
                    "description": "bla",
                    "name": "bla3",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": [
                    "trans1a",
                    "trans1b",
                    "trans1c",
                    "trans2",
                    "trans3",
                ],
                "gene_stable_id": ["fake1", "fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "1", "2"],
                "strand": [1, 1, 1, -1, -1],
                "start": [3100, 3000, 4850, 4910, 4900],
                "stop": [4900, 4000, 4950, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4850, 4950)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    g = genome.genes["fake1"]
    one = g.exons_merged
    assert (one[0] == [3000]).all()
    assert (one[1] == [4950]).all()
    g = genome.genes["fake2"]
    two = g.exons_merged
    assert (two[0] == [4910, 5100]).all()
    assert (two[1] == [5000, 5400]).all()
    four = genome.genes["fake4"].exons_merged
    assert len(four[0]) == 0


def test_get_gene_exons_protein_coding_merged():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                    "biotype": "protein_coding",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                    "biotype": "protein_coding",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                    "biotype": "lincRNA",
                },
                {
                    "stable_id": "fake4",
                    "chr": "3",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                    "biotype": "noncoding_rna",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": [
                    "trans1a",
                    "trans1b",
                    "trans2",
                    "trans3",
                    "trans4",
                ],
                "biotype": [
                    "protein_coding",
                    "whatever",
                    "protein_coding",
                    "lincRNA",
                    "non_coding_rna",
                ],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3", "fake4"],
                "chr": ["1", "1", "1", "2", "3"],
                "strand": [1, 1, -1, -1, -1],
                "start": [3100, 3000, 4910, 4900, 4900],
                "stop": [3200, 4000, 5400, 5400, 5400],
                "exons": [
                    [(3100, 3200)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                    [],
                ],
            }
        ),
    )
    g = genome.genes["fake1"]
    one = g.exons_protein_coding_merged
    assert (one[0] == [3100]).all()
    assert (one[1] == [3200]).all()

    g = genome.genes["fake2"]
    two = g.exons_protein_coding_merged
    assert (two[0] == [4910, 5100]).all()
    assert (two[1] == [5000, 5400]).all()

    g = genome.genes["fake3"]
    three = g.exons_protein_coding_merged
    assert (three[0] == [4900]).all()
    assert (three[1] == [5400]).all()
    four = genome.genes["fake4"].exons_protein_coding_merged
    assert len(four[0]) == 0


def test_get_gene_exons_protein_coding_overlapping():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                    "biotype": "protein_coding",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                    "biotype": "protein_coding",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                    "biotype": "lincRNA",
                },
                {
                    "stable_id": "fake4",
                    "chr": "3",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                    "biotype": "noncoding_rna",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": [
                    "trans1a",
                    "trans1b",
                    "trans2a",
                    "trans2b",
                    "trans3",
                    "trans4",
                ],
                "biotype": [
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "non_coding",
                    "lincRNA",
                    "non_coding_rna",
                ],
                "gene_stable_id": [
                    "fake1",
                    "fake1",
                    "fake2",
                    "fake2",
                    "fake3",
                    "fake4",
                ],
                "chr": ["1", "1", "1", "2", "2", "3"],
                "strand": [1, 1, -1, -1, -1, -1],
                "start": [3100, 3000, 4910, 4950, 4900, 4900],
                "stop": [3200, 4000, 5400, 5300, 5400, 5400],
                "exons": [
                    [(3100, 3200)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4950, 5300)],
                    [(4900, 5400)],
                    [],
                ],
            }
        ),
    )
    g = genome.genes["fake1"]
    one = g.exons_protein_coding_overlapping
    print(one)
    assert (one[0] == [3000, 3100, 3300, 3750]).all()
    assert (one[1] == [3500, 3200, 3330, 4000]).all()

    g = genome.genes["fake2"]
    two = g.exons_protein_coding_overlapping
    assert (two[0] == [4910, 5100]).all()
    assert (two[1] == [5000, 5400]).all()

    g = genome.genes["fake3"]
    three = g.exons_protein_coding_overlapping
    assert (three[0] == [4900]).all()
    assert (three[1] == [5400]).all()
    four = genome.genes["fake4"].exons_protein_coding_overlapping
    assert len(four[0]) == 0


def test_gene_exons_overlapping():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
                {
                    "stable_id": "fake4",
                    "chr": "2",
                    "strand": -1,
                    "tss": 6400,
                    "tes": 5900,
                    "description": "bla",
                    "name": "bla3",
                },
            ]
        ),
        # {transcript_stable_id, gene_stable_id, strand, start, end, exons},
        df_transcripts=pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4910, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    one = genome.genes["fake1"].exons_overlapping
    print(one)
    assert (one[0] == [3000, 3100, 3300, 3750]).all()
    assert (one[1] == [3500, 4900, 3330, 4000]).all()
    assert len(genome.genes["fake4"].exons_overlapping[0]) == 0


def test_gene_tss_tes():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
                {
                    "stable_id": "fake4",
                    "chr": "2",
                    "strand": -1,
                    "tss": 6400,
                    "tes": 5900,
                    "description": "bla",
                    "name": "bla3",
                },
            ]
        ),
        pd.DataFrame(
            {
                "transcript_stable_id": ["trans1a", "trans1b", "trans2", "trans3"],
                "gene_stable_id": ["fake1", "fake1", "fake2", "fake3"],
                "chr": ["1", "1", "1", "2"],
                "strand": [1, 1, -1, -1],
                "start": [3100, 3000, 4910, 4900],
                "stop": [4900, 4000, 5400, 5400],
                "exons": [
                    [(3100, 4900)],
                    [(3000, 3500), (3300, 3330), (3750, 4000)],
                    [(4910, 5000), (5100, 5400)],
                    [(4900, 5400)],
                ],
            }
        ),
    )
    assert genome.genes["fake1"].tss == 3000
    assert genome.genes["fake1"].tes == 4900
    assert genome.genes["fake2"].tss == 5400
    assert genome.genes["fake2"].tes == 4900


def test_name_to_gene_id():
    genome = DummyGenome(
        pd.DataFrame(
            [
                {
                    "gene_stable_id": "fake1",
                    "chr": "1",
                    "strand": 1,
                    "tss": 3000,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla1",
                },
                {
                    "gene_stable_id": "fake2",
                    "chr": "1",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla2",
                },
                {
                    "gene_stable_id": "fake3",
                    "chr": "2",
                    "strand": -1,
                    "tss": 5400,
                    "tes": 4900,
                    "description": "bla",
                    "name": "bla3",
                },
                {
                    "gene_stable_id": "fake4",
                    "chr": "2",
                    "strand": -1,
                    "tss": 6400,
                    "tes": 5900,
                    "description": "bla",
                    "name": "bla4",
                },
            ]
        )
    )
    assert genome.name_to_gene_ids("bla1") == set(["fake1"])
    assert genome.name_to_gene_ids("bla2") == set(["fake2"])
    assert genome.name_to_gene_ids("bla3") == set(["fake3"])
    assert genome.name_to_gene_ids("bla4") == set(["fake4"])


def test_get_reads_in_exon():
    import mbf_sampledata
    import pysam

    genome = mbf_sampledata.get_human_22_fake_genome()
    bam = pysam.Samfile(
        mbf_sampledata.get_sample_path("mbf_align/rnaseq_spliced_chr22.bam")
    )
    g = genome.genes["ENSG00000128228"]
    reads = g.get_reads_in_exons(bam)
    assert reads
    start = 21642302 - 1
    stop = 21644299
    for r in reads:
        ov = r.get_overlap(start, stop)
        assert ov > 0


def test_transcript_coordinate_translation_fwd():
    import mbf_sampledata

    genome = mbf_sampledata.get_human_22_fake_genome()
    tr = genome.transcripts["ENST00000455558"]
    map = tr.coordinate_translations
    assert len(map) == 657
    # first exon
    for ii in range(0, 47):
        assert (25565107 - 1 + ii) in map

    reference = [
        (25565107, 25565153, 47),
        (25604377, 25604453, 77),
        (25661576, 25661677, 102),
        (25663630, 25663704, 75),
        (25667739, 25667800, 62),
        (25672296, 25672347, 52),
        (25674437, 25674528, 92),
        (25678816, 25678915, 100),
        (25684519, 25684568, 50),
    ]
    offset = 0
    should = []
    for genome_start, genome_stop, length in reference:
        for ii in range(offset, offset + length):
            should.append(genome_start - 1 + ii - offset)
        offset += length
    ok = map == should
    assert ok  # pytest spends an eternity formating the diff otherwise on failure..

    assert map[-1] == (25684568 - 1)


def test_transcript_coordinate_translation_reverse():
    import mbf_sampledata

    genome = mbf_sampledata.get_human_22_fake_genome()
    tr = genome.transcripts["ENST00000420242"]
    map = tr.coordinate_translations
    assert len(map) == 525

    print(tr.exons)
    assert map[0] == (26512496 - 1)
    assert map[102] == (26512394 - 1)
    assert map[103] == (26512175 - 1)
    for ii in range(77):
        assert map[103 + ii] == (26512175 - 1 - ii)
    assert map[103 + 77] == (26510692 - 1)
    assert map[103 + 77 + 76 - 1] == (26510617 - 1)
    assert map[-1] == (26506878 - 1)
    assert len(set(map)) == len(map)
