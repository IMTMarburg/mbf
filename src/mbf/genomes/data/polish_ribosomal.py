import pandas as pd
import csv
import re
import collections
from pathlib import Path


def fix_chromosome(c):
    c = c[3:]
    if c == "M":
        return "MT"
    if len(c) <= 2:
        return c
    if "gl" in c:
        return c[c.find("gl") :].upper() + ".1"
    return "remove"


def parse_attrs(attrs):
    attrs = re.findall('([^ ]+) "([^"]+)";', attrs)
    return dict(attrs)


def format_attrs(attrs):
    return " ".join(['%s "%s"; ' % x for x in attrs.items()])


for input_filename in Path(__file__).parent.glob("ribosomal_genes_*.gtf.gz"):
    if "full" in input_filename.name:
        continue
    print("parsing", input_filename)
    df = pd.read_csv(input_filename, header=None, sep="\t", comment="#")
    df.columns = [
        "chr",
        "source",  # constant
        "kind",  # constant
        "start",  # uscs database data is 0 based...
        "stop",
        "score",
        "strand",
        "ignored",
        "annotation",
    ]
    df["chr"] = [fix_chromosome(x) for x in df["chr"]]
    df = df[df['chr'] != 'remove']

    if (df.kind != "exon").any():
        raise ValueError("non exon found - unexpected")
    # ok, this file is all exons, we need to create one gene, one transcript per row basicially,
    # no splicing here,
    result = []
    dedup = collections.Counter()
    for idx, row in df.iterrows():
        attrs = parse_attrs(row["annotation"])
        gene_id = attrs["gene_id"] + "_" + str(dedup[attrs["gene_id"]])
        dedup[attrs["gene_id"]] += 1
        gene_name = gene_id
        biotype = "tRNA" if "tRNA" in gene_id else "rRNA"
        attrs = {
            "gene_id": gene_id,
            "gene_name": gene_name,
            "gene_source": "ucsc_rmsk",
            "gene_biotype": biotype,
        }
        gene = {
            "chr": row["chr"],
            "source": "ucsc_rmsk",
            "kind": "gene",
            "start": row["start"],
            "stop": row["stop"],
            "score": row["score"],
            "strand": row["strand"],
            "ignored": ".",
            "annotation": format_attrs(attrs),
        }
        result.append(gene)

        attrs.update(
            {
                "transcript_id": gene_id + "_tr",
                "transcript_name": gene_id,
                "transcript_source": "ucsc_rmsk",
                "transcript_biotype": biotype,
            }
        )
        transcript = {
            "chr": row["chr"],
            "source": "ucsc_rmsk",
            "kind": "transcript",
            "start": row["start"],
            "stop": row["stop"],
            "score": row["score"],
            "strand": row["strand"],
            "ignored": ".",
            "annotation": format_attrs(attrs),
        }
        result.append(transcript)

        attrs.update({"exon_id": gene_id + "_exon"})
        exon = {
            "chr": row["chr"],
            "source": "ucsc_rmsk",
            "kind": "exon",
            "start": row["start"],
            "stop": row["stop"],
            "score": row["score"],
            "strand": row["strand"],
            "ignored": ".",
            "annotation": format_attrs(attrs),
        }
        result.append(exon)
    output_df = pd.DataFrame(result)[df.columns]
    output_df.to_csv(
        input_filename.with_name(input_filename.name + ".full.gtf.gz"),
        sep="\t",
        header=False,
        quoting=csv.QUOTE_NONE,
        index=False,
    )
