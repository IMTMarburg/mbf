from ..annotator import Annotator
import pandas as pd
import pypipegraph as ppg
import numpy
import collections
from pathlib import Path
import dppd
import dppd_plotnine  # noqa: F401

dp, X = dppd.dppd()

_GenomicPosition_Stiewe_GenomicDistribution = {}


class GenomicPositionStiewe(Annotator):
    """This annotator adds
    promotor (+-5 kb around TSS), genebody (rest of TSS-TES) or distal (either -25kb to -5 kb or TES to TES+25kb), intergenic
    #######'close to TSS', 'intron', 'exon',
    There is only one case, TSS trumps everything, than exon, then intron, then upstream, then intragenic
    """

    _plot_jobs = {}

    def __init__(self, tss_distance=5000, distal_distance=25000):
        self.columns = [
            "GenomicPosition TSS=%i bp, distal=%i bp" % (tss_distance, distal_distance)
        ]
        self.tss_distance = tss_distance
        self.distal_distance = distal_distance
        #        self.upstream_distance = upstream_distance
        self.column_properties = {
            self.columns[0]: {
                "description": "Adds one of 'promotor', 'gene_body', 'TSS distal', 'TES distal', 'intergenic', earlier entries beat later entries"
            }
        }
        Annotator.__init__(self)
        pass

    def deps(self, genomic_regions):
        from ..genes import Genes
        from ..transcripts import Transcripts

        all_genes = Genes(genomic_regions.genome)
        all_transcripts = Transcripts(genomic_regions.genome)
        tss = all_transcripts.regions_tss()
        tes = all_transcripts.regions_tes()
        introns = all_genes.regions_introns()
        exons = all_genes.regions_exons_merged()
        return [
            all_genes.load(),
            tss.load(),
            tes.load(),
            introns.load(),
            exons.load(),
            # don't need a ParameterInvariant, since it will change the column & get_name
            # ppg.ParameterInvariant(genomic_regions.name + self.columns, (self.tss_distance, ))
        ]

    def calc_ddf(self, genomic_regions):
        result = []
        from ..transcripts import Transcripts
        from ..genes import Genes

        all_transcripts = Transcripts(genomic_regions.genome)
        tss = all_transcripts.regions_tss()
        tes = all_transcripts.regions_tes()
        genes = Genes(genomic_regions.genome)
        tss_distance = self.tss_distance  # the distance we consider close to a tss/tes
        # upstream_distance = self.upstream_distance
        for dummy_idx, row in genomic_regions.df[["chr", "start", "stop"]].iterrows():
            here = "intergenic"
            overlaps_tss = tss.has_overlapping(
                row["chr"], max(0, row["start"] - tss_distance), row["stop"] + tss_distance
            )
            if overlaps_tss:
                here = "promotor"
            else:
                overlap_gene_body = genes.has_overlapping(
                    row["chr"], row["start"], row["stop"]
                )
                if overlap_gene_body:
                    here = "gene body"
                else:
                    overlap_distal = tss.has_overlapping(
                        row["chr"],
                        max(0, row["start"] - self.distal_distance),
                        row["stop"] + self.distal_distance,
                    )
                    if overlap_distal:
                        here = "TSS distal"
                    else:
                        overlap_upstream_distal = tes.has_overlapping(
                            row["chr"],
                            max(0, row["start"] - self.distal_distance),
                            row["stop"] + self.distal_distance,
                        )
                        if overlap_upstream_distal:
                            here = "TES distal"
            result.append(here)
        return pd.DataFrame({self.columns[0]: result})

    def plot(self, genomic_regions):
        of = genomic_regions.result_dir / (
            "GenomicPosition_Stiewe_%i_%i.pdf"
            % (self.tss_distance, self.distal_distance)
        )
        if not of in self._plot_jobs:
            of.parent.mkdir(exist_ok=True, parents=True)

            def calc():
                counts = collections.Counter()
                for entry in genomic_regions.df[self.columns[0]].values:
                    counts[entry] += 1
                return pd.DataFrame({"Region": counts.keys(), "Count": counts.values()})

            def plot(df):
                if len(df) == 0:
                    return (
                        dp(pd.DataFrame({"text": "no regions"}))
                        .p9()
                        .add_text(_x=0, y_=0, label="text")
                    ).pd
                df["Source"] = [genomic_regions.name] * len(df)
                p = (
                    dp(df)
                    .p9()
                    .add_bar("Source", "Count", fill="Region", position="stack")
                ).pd
                return p

            job = ppg.PlotJob(of, calc, plot).depends_on(
                genomic_regions.add_annotator(self)
            )
            ppg.Job.depends_on(
                job, self.calc_genomic_distribution(genomic_regions.genome)
            )
            self._plot_jobs[of] = job
        return self._plot_jobs[of]

    def calc_genomic_distribution(self, genome):
        from ..genes import Genes
        from ..transcripts import Transcripts

        all_genes = Genes(genome)
        all_transcripts = Transcripts(genome)
        # tes = all_transcripts.regions_tes()
        tss = all_transcripts.regions_tss()
        introns = all_genes.regions_introns()
        # transcripts = all_genes.regions_body()
        exons = all_genes.regions_exons_merged()

        def calc():
            counts = [
                "promotor",
                "gene body",
                "TES distal",
                "TSS distal",
                "intergenic",
            ]
            counts = dict((k, 0) for k in counts)
            total_length = 0
            for chr, chr_length in genome.get_chromosome_lengths().items():
                total_length += chr_length
                array = numpy.zeros((chr_length,), dtype=numpy.uint8)
                code_tss = 1
                code_tes_distal = 2
                code_tss_distal = 3
                code_gene_body = 4
                code_intergenic = 0
                for _idx, (start, stop, strand) in all_transcripts.df[
                    all_transcripts.df["chr"] == chr
                ][["start", "stop", "strand"]].iterrows():
                    start = int(start)
                    stop = int(stop)
                    if strand == 1:
                        array[start - self.distal_distance : stop] = code_tss_distal
                        array[
                            start : stop + self.distal_distance
                        ] = code_tes_distal
                    else:
                        array[
                            start - self.distal_distance : stop
                        ] = code_tes_distal
                        array[start : stop + self.distal_distance] = code_tss_distal
                for _idx, (start, stop) in all_transcripts.df[
                    all_transcripts.df["chr"] == chr
                ][["start", "stop"]].iterrows():
                    start = int(start)
                    stop = int(stop)
                    array[start:stop] = code_gene_body
                for _idx, (start, stop) in tss.df[tss.df["chr"] == chr][
                    ["start", "stop"]
                ].iterrows():
                    start = int(start)
                    stop = int(stop)
                    array[
                        start - self.tss_distance : stop + self.tss_distance
                    ] = code_tss
                counts["intergenic"] += numpy.sum(array == code_intergenic)
                counts["TES distal"] += numpy.sum(array == code_tes_distal)
                counts["TSS distal"] += numpy.sum(array == code_tss_distal)
                counts["gene body"] += numpy.sum(array == code_gene_body)
                counts["promotor"] += numpy.sum(array == code_tss)

            if sum(counts.values()) != total_length:
                import pprint

                pprint.pprint(counts)
                raise ValueError(
                    "Sum of counts was not correct was: %i, but we have %i bp"
                    % (sum(counts.values()), total_length)
                )
            return pd.DataFrame(
                {"Region": list(counts.keys()), "Count": list(counts.values())}
            )

        def store(value):
            _GenomicPosition_Stiewe_GenomicDistribution[
                genome.name, self.tss_distance, self.distal_distance
            ] = value

        fn = Path(
            "cache/GenomicRegions/Annotator_GenomicPosition_Stiewe/%s_%i_%i.dat"
            % (genome.name, self.tss_distance, self.distal_distance)
        )
        fn.parent.mkdir(exist_ok=True, parents=True)
        return ppg.CachedDataLoadingJob(fn, calc, store).depends_on(
            [
                all_genes.load(),
                tss.load(),
                introns.load(),
                exons.load(),
                all_transcripts.load(),
            ]
        )

    def plot_bars_vs_genomic_distribution(self, gr):
        of = gr.result_dir / (
            "GenomicPosition_Stiewe_Global_vs_local_bar_%i_%ibp.pdf"
            % (self.tss_distance, self.distal_distance)
        )

        def calc():
            data_genomic = _GenomicPosition_Stiewe_GenomicDistribution[
                gr.genome.name, self.tss_distance, self.distal_distance
            ]
            data_local = self.plot(gr).data_.copy()
            local_total = float(numpy.sum(data_local["Count"]))
            genomic_total = float(numpy.sum(data_genomic["Count"]))
            data_local = data_local.assign(
                Relative=[(x / local_total) for x in data_local["Count"]]
            )
            data_genomic = data_genomic.assign(
                Relative=[(x / genomic_total) for x in data_genomic["Count"]]
            )

            data_local = data_local.assign(Source=[gr.name] * len(data_local))
            data_genomic = data_genomic.assign(Source=["Genomic"] * len(data_genomic))
            # if (data_genomic == data_local).all().all(): # how could that even happen... they have Source being different
            #     raise ValueError("identical")
            return pd.concat([data_local, data_genomic], axis=0)

        def plot(df):
            # for ii, val in enumerate(
            #     ["intergenic", "distal", "gene body", "promotor", "upstream distal"]
            # ):
            #     # TODO: pandas check
            #     df[df["Region"] == val, "Region"] = "%i. %s" % (ii, val)
            # df = df.sort_by("Region")
            p = (
                dp(df)
                #.categorize( 'Region', ['intergenic', 'intron', 'exon', 'tss', 'upstream'])
                .p9()
                .theme(legend_position='bottom')
                .add_bar("Source", "Relative", position="stack", fill="Region")
                .coord_flip()
                #.scale_fill_brewer(palette="Set3")
            ).pd
            return p

        of.parent.mkdir(exist_ok=True, parents=True)
        return (
            ppg.PlotJob(of, calc, plot)
            .depends_on(self.plot(gr).cache_job)
            .depends_on(self.calc_genomic_distribution(gr.genome))
        )
