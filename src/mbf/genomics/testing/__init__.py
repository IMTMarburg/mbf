from mbf_genomes import HardCodedGenome


def MockGenome(  # noqa: C901
    df_genes, df_transcripts=None, chr_lengths=None, df_genes_meta=None
):  # noqa: C901
    if chr_lengths is None:
        chr_lengths = {
            "1": 100_000,
            "2": 200_000,
            "3": 300_000,
            "4": 400_000,
            "5": 500_000,
        }

    df_genes = df_genes.rename(columns={"stable_id": "gene_stable_id"})
    if not "start" in df_genes.columns:
        starts = []
        stops = []
        if not "strand" in df_genes:
            tes_larger = df_genes["tes"] > df_genes["tss"]
            strand = tes_larger.replace({True: 1, False: -1})
            df_genes = df_genes.assign(strand=strand)
        for idx, row in df_genes.iterrows():
            starts.append(min(row["tss"], row["tes"]))
            stops.append(max(row["tss"], row["tes"]))
        df_genes = df_genes.assign(start=starts, stop=stops)
    if not "biotype" in df_genes.columns:
        df_genes = df_genes.assign(biotype="protein_coding")
    if not "name" in df_genes.columns:
        df_genes = df_genes.assign(name=df_genes.gene_stable_id)
    df_genes = df_genes.sort_values(["chr", "start"])
    df_genes = df_genes.set_index("gene_stable_id")
    if not df_genes.index.is_unique:  # pragma: no cover
        raise ValueError("gene_stable_ids not unique")
    if df_transcripts is not None and len(df_transcripts):
        if not "transcript_stable_id" in df_transcripts.columns:
            print(df_transcripts.columns)
            df_transcripts = df_transcripts.assign(
                transcript_stable_id=df_transcripts["name"]
            )
        if not "biotype" in df_transcripts.columns:
            df_transcripts = df_transcripts.assign(biotype="protein_coding")
        if not "name" in df_transcripts.columns:
            df_transcripts = df_transcripts.assign(
                name=df_transcripts.transcript_stable_id
            )
        if "exons" in df_transcripts.columns:
            if len(df_transcripts["exons"].iloc[0]) == 3:  # pragma: no cover
                df_transcripts = df_transcripts.assign(
                    exons=[(x[0], x[1]) for x in df_transcripts["exons"]]
                )
            df_transcripts = df_transcripts.assign(
                exon_stable_ids=[
                    "exon_%s_%i" % (idx, ii)
                    for (ii, idx) in enumerate(df_transcripts["exons"])
                ]
            )
        stops = []
        if not "strand" in df_transcripts:  # pragma: no cover
            df_transcripts = df_transcripts.assign(strand=1)
        if not "tss" in df_transcripts:  # pragma: no branch
            tss = []
            tes = []
            for tr, row in df_transcripts.iterrows():
                if row["strand"] == 1:
                    tss.append(min((x[0] for x in row["exons"])))
                    tes.append(max((x[1] for x in row["exons"])))
                else:
                    tss.append(max((x[1] for x in row["exons"])))
                    tes.append(min((x[0] for x in row["exons"])))
            df_transcripts = df_transcripts.assign(tss=tss, tes=tes)
        if not "start" in df_transcripts:
            starts = []
            stops = []
            for idx, row in df_transcripts.iterrows():
                starts.append(min(row["tss"], row["tes"]))
                stops.append(max(row["tss"], row["tes"]))
            df_transcripts = df_transcripts.assign(start=starts, stop=stops)
        df_transcripts = df_transcripts.set_index("transcript_stable_id")
        if not df_transcripts.index.is_unique:  # pragma: no cover
            raise ValueError("transcript_stable_ids not unique")
    result = HardCodedGenome("dummy", chr_lengths, df_genes, df_transcripts, None)
    result.sanity_check_genes(df_genes)
    if df_transcripts is not None and len(df_transcripts):
        result.sanity_check_transcripts(df_transcripts)
    if df_genes_meta is not None:
        result.df_genes_meta = df_genes_meta
    return result
