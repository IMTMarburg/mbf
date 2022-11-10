import pandas as pd
import subprocess
import io
import dppd
from pathlib import Path
import mbf.functional
from statsmodels.stats.multitest import multipletests

dp, X = dppd.dppd()


def get_gene_name(genome, x):
    names = genome.gene_id_to_name(x)
    if not names:
        return x
    return names[0]


def multi_hypergeom_test(
    genome, query_set, function_gene_groups_or_list_of_such=None, background_set=None
):
    """Test a query set against multiple sets from functional.databases.
    Returns a pandas.DataFrame(group, set, benjamini, p, overlap_count,
    sorted by benjamini

    """
    if function_gene_groups_or_list_of_such is None:
        from . import databases

        function_gene_groups_or_list_of_such = databases.get_default_groups()

    with open("input.gmt", mode="w") as handle:
        for func_group in function_gene_groups_or_list_of_such:
            sets = func_group.get_sets(genome)
            for geneset_name, genes in sets.items():
                row = "%s:::%s\tNA\t" % (func_group.name, geneset_name)
                handle.write(row)
                try:
                    handle.write("\t".join(set(genes)))
                    handle.write("\n")
                except TypeError:
                    print(func_group.name)
                    print(geneset_name)
                    print(genes)
        handle.flush()

        # that's from the FunctionalHypergeomTest repo
        # and does functional testing at the speed of Rust
        p = subprocess.Popen(
            [
                "functional_hypergeom_test",
                "test",
                "--reference_sets",
                handle.name,
                "--full",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
        )
        stdout, stderr = p.communicate(("\n".join(query_set)).encode("utf-8"))
        if p.returncode != 0:
            raise ValueError(
                "functional_hypergeom_test error return: %i %s %s"
                % (p.returncode, stdout, stderr)
            )
        df = pd.read_csv(io.StringIO(stdout.decode("utf-8")), sep="\t", index_col=False)
        if len(df):
            df = (
                dp(df).seperate("Set_name", ["group", "set"], sep=":::", remove=True)
            ).pd
        df = dp(df).rename(columns={"P_value": "p", "Overlap": "overlap count"}).pd
        if len(df):
            fdr = multipletests(df["p"], method="fdr_bh")[1]
        else:
            fdr = []
        df = df.assign(benjamini=fdr)
        df = df.sort_values("benjamini")
        return df
