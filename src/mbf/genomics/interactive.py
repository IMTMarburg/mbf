"""interactive helpers"""
import numpy as np
import pandas as pd
import dppd
dp, X = dppd.dppd()


def read_default_genes(
    fn="results/Genes/protein_coding/protein_coding.tsv", kind="auto"
):
    """I often have to read the genes df to CPM/TPM,
    and it's annoying"""
    import pandas as pd
    import dppd

    dp, X = dppd.dppd()
    df = pd.read_csv(fn, sep="\t")
    exon_columns = [
        x
        for x in df.columns
        if x.startswith("Exon, protein coding, stranded smart tag count ")
    ]
    if kind == "auto":
        cpm = any((x.endswith("CPM") for x in exon_columns))
        tpm = any((x.endswith("TPM") for x in exon_columns))
        if cpm and tpm:
            raise ValueError("both TPM and CPM present, pass kind = 'CPM' or 'TPM'")
        if cpm:
            kind = "CPM"
        elif tpm:
            kind = "TPM"
        else:
            raise ValueError("No exon columns CPM/TPM found")
    if kind == "raw":
        xpm_columns = [x for x in exon_columns if x.endswith("_STAR")]
    else:
        xpm_columns = [x for x in exon_columns if x.endswith(kind)]
    out = (
        dp(df)
        .set_index(["gene_stable_id", "name"])[xpm_columns]
        .reset_columns(
            lambda x: x.replace("Exon, protein coding, stranded smart tag count ", "")
            .replace(kind, "")
            .replace("_STAR", "")
            .strip()
        )
    ).pd
    return out


def ma(df, a, b, epsilon=0.1):
    """quick MA plot wrapper"""
    la = np.log2(df[a] + epsilon)
    lb = np.log2(df[b] + epsilon)
    A = (la + lb) / 2
    M = la - lb
    pdf = pd.DataFrame({"A": A, "M": M}, index=df.index)
    return (
        dp(pdf)
        .p9()
        .add_scatter("A", "M")
        .add_hline(0, _color="blue")
        .add_hline(1, _color="blue")
        .add_hline(-1, _color="blue")
    ).pd


__all__ = ['ma','read_default_genes']
