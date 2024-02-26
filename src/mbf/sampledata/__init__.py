import pandas as pd
from pathlib import Path

__version__ = "0.2"


def get_sample_data(fn) -> str:
    here = Path(__file__).parent
    return str((here / "data" / fn).absolute())


def get_sample_path(fn) -> Path:
    here = Path(__file__).parent
    return (here / "data" / fn).absolute()


def get_human_22_fake_genome():
    from mbf.genomics.testing import MockGenome
    import gzip

    genes = pd.read_parquet(
        gzip.GzipFile(get_sample_data(Path("mbf_align/hs_22_genes.parquet.gz")))
    ).reset_index()
    tr = pd.read_parquet(
        gzip.GzipFile(get_sample_data(Path("mbf_align/hs_22_transcripts.parquet.gz")))
    ).reset_index()
    return MockGenome(df_genes=genes, df_transcripts=tr, chr_lengths={"22": 50_818_468})


def get_Candidatus_carsonella_ruddii_pv(name=None, **kwargs):
    """A FilebasedGenome used by other libraries for their tests"""
    from mbf.genomes import FileBasedGenome

    if name is None:  # pragma: no cover
        name = "Candidatus_carsonella"
    return FileBasedGenome(
        name,
        get_sample_path(
            "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.dna.toplevel.fa.gz"
        ),
        get_sample_path(
            "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.42.gtf.gz"
        ),
        get_sample_path(
            "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.cdna.all.fa.gz"
        ),
        get_sample_path(
            "mbf_genomes/Candidatus_carsonella_ruddii_pv.ASM1036v1.pep.all.fa.gz"
        ),
        **kwargs,
    )


def get_pasilla_data_subset():
    from mbf.genomics import DelayedDataFrame
    import numpy as np

    pasilla_data = pd.read_csv(
        get_sample_path("mbf_comparisons/pasillaCount_deseq2.tsv.gz"), sep=" "
    )
    # pasilla_data = pasilla_data.set_index('Gene')
    pasilla_data.columns = [str(x) for x in pasilla_data.columns]
    treated = [x for x in pasilla_data.columns if x.startswith("treated")]
    untreated = [x for x in pasilla_data.columns if x.startswith("untreated")]
    pasilla_data = pasilla_data.assign(
        abslog2FC=np.log2(pasilla_data[treated].mean(axis=1) + 0.1)
        - np.log2(pasilla_data[untreated].mean(axis=1) + 0.1)
    )

    pasilla_data = pasilla_data[pasilla_data[treated + untreated].mean(axis=1) > 50]
    pasilla_data = pasilla_data.sort_values("abslog2FC", ascending=False)[:1000]
    pasilla_data = DelayedDataFrame("pasilla", pasilla_data)
    return pasilla_data, treated, untreated
