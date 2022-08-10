from ..regions import GenomicRegions

import numpy as np
import pypipegraph as ppg
import pandas as pd
from pathlib import Path
from mbf_externals.util import lazy_method

_exon_regions_overlapping_cache = {}
_genes_per_genome_singletons = {}


class Genes(GenomicRegions):
    def __new__(cls, genome, alternative_load_func=None, *args, **kwargs):
        """Make sure that Genes for a full genome (ie. before filtering) are singletonic. Ie.
        you can always safely call Genes(my_genome), without worrying about duplicate objects"""
        if alternative_load_func is None:
            if ppg.util.global_pipegraph:
                if not hasattr(
                    ppg.util.global_pipegraph, "_genes_per_genome_singletons"
                ):
                    ppg.util.global_pipegraph._genes_per_genome_singletons = {}
                singleton_dict = ppg.util.global_pipegraph._genes_per_genome_singletons
            else:
                singleton_dict = _genes_per_genome_singletons
            if not genome in singleton_dict:
                singleton_dict[genome] = GenomicRegions.__new__(cls)
            return singleton_dict[genome]
        else:
            return GenomicRegions.__new__(cls)

    def __init__(
        self,
        genome,
        alternative_load_func=None,
        name=None,
        dependencies=None,
        result_dir=None,
        sheet_name=None,
        vid=None,
    ):
        if hasattr(self, "_already_inited"):
            if (
                alternative_load_func is not None
                and alternative_load_func != self.genes_load_func
            ):  # pragma: no cover -
                # this can only happen iff somebody starts to singletonize the Genes
                # with loading_functions - otherwise the duplicate object
                # checker will kick in first
                raise ValueError(
                    "Trying to define Genes(%s) twice with different loading funcs"
                    % self.name
                )
            pass
        else:
            if name is None:
                if alternative_load_func:
                    raise ValueError(
                        "If you pass in an alternative_load_func you also need to specify a name"
                    )
                name = "Genes_%s" % genome.name
                self.top_level = True
            else:
                self.top_level = False
            if alternative_load_func is None:
                load_func = lambda: genome.df_genes.reset_index()  # noqa: E731
            else:
                load_func = alternative_load_func

            self.genes_load_func = load_func
            if result_dir:
                pass
            elif sheet_name:
                result_dir = result_dir or Path("results") / "Genes" / sheet_name / name
            else:
                result_dir = result_dir or Path("results") / "Genes" / name
            # result_dir = Path(result_dir).absolute() # must not be absolute in ppg2

            self.column_properties = {
                "chr": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "On which chromosome (or contig) the gene is loacted",
                },
                "start": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "Left most position of this gene",
                },
                "stop": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "Right most position of this gene",
                },
                "tss": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "Position of Transcription start site on the chromosome",
                },
                "tes": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "Position of the annotated end of transcription",
                },
                "gene_stable_id": {
                    "user_visible": True,
                    "priority": -1000,
                    "index": True,
                    "description": "Unique identification name for this gene",
                },
                "name": {
                    "user_visible": True,
                    "priority": -999,
                    "index": True,
                    "description": "Offical name for this gene",
                    "nocase": True,
                },
                "strand": {
                    "user_visible": True,
                    "priority": -998,
                    "description": "Which is the coding strand (1 forward, -1 reverse)",
                },
                "description": {
                    "user_visible": True,
                    "priority": -998.5,
                    "description": "A short description of the gene",
                },
            }
            GenomicRegions.__init__(
                self,
                name,
                self._load,
                dependencies if dependencies is not None else [],
                genome,
                on_overlap="ignore",
                result_dir=result_dir,
                sheet_name=sheet_name,
                summit_annotator=False,
            )
            if self.load_strategy.build_deps:
                deps = [genome.download_genome()]
                # function invariant for _load is done by GR.__init__
                deps.append(genome.job_genes())
                deps.append(genome.job_transcripts())
                self.load_strategy.load().depends_on(deps)
            self._already_inited = True
            self.vid = vid

    def __str__(self):
        return "Genes(%s)" % self.name

    def __repr__(self):
        return "Genes(%s)" % self.name

    def register(self):
        pass

    def _load(self):
        """Load func
        """
        if hasattr(self, "df"):  # pragma: no cover
            return
        df = self.genes_load_func()

        if not isinstance(df, pd.DataFrame):
            raise ValueError(
                "GenomicRegion.loading_function must return a DataFrame, was: %s"
                % type(df)
            )
        for col in self.get_default_columns():
            if not col in df.columns:
                func_filename = self.genes_load_func.__code__.co_filename
                func_line_no = self.genes_load_func.__code__.co_firstlineno
                raise ValueError(
                    "%s not in dataframe returned by GenomicRegion.loading_function %s %s - it did return %s"
                    % (col, func_filename, func_line_no, df.columns)
                )
        allowed_chromosomes = set(self.genome.get_chromosome_lengths().keys())
        if len(df):
            for chr in df["chr"]:
                if not chr in allowed_chromosomes:
                    raise ValueError(
                        "Invalid chromosome found when loading %s: %s, expected one of: %s\nLoading func was %s"
                        % (
                            self.name,
                            repr(chr),
                            sorted(allowed_chromosomes),
                            self.genes_load_func,
                        )
                    )
            if not np.issubdtype(df["tss"].dtype, np.integer):
                raise ValueError(
                    "tss needs to be an integer, was: %s" % df["tss"].dtype
                )
            if not np.issubdtype(df["tes"].dtype, np.integer):
                raise ValueError(
                    "tes needs to be an integer, was: %s" % df["tes"].dtype
                )
            # df = self.handle_overlap(df) Genes don't care about overlap
            if not "start" in df.columns:
                df = df.assign(start=np.min([df["tss"], df["tes"]], 0))
            if not "stop" in df.columns:
                df = df.assign(stop=np.max([df["tss"], df["tes"]], 0))
        else:
            df = df.assign(
                start=np.array([], dtype=np.int32), stop=np.array([], dtype=np.int32)
            )

        if (df["start"] > df["stop"]).any():
            raise ValueError(
                "Genes.loading_function returned a negative interval:\n %s"
                % df[df["start"] > df["stop"]].head()
            )
        self.df = df.sort_values(["chr", "start"], ascending=[True, True]).reset_index(
            drop=True
        )  # since we don't call handle_overlap
        # enforce column order
        cols = ["gene_stable_id", "chr", "start", "stop", "strand", "tss", "tes"]
        for x in df.columns:
            if not x in cols:
                cols.append(x)
        df = df[cols]
        return df

    def get_default_columns(self):
        return ("chr", "tss", "tes", "gene_stable_id")

    def _new_for_filtering(self, new_name, load_func, dependencies, **kwargs):
        """When filtering, a new object of this class is created.
        To pass it the right options from the parent, overwrite this
        """
        return Genes(self.genome, load_func, new_name, dependencies, **kwargs)

    def overlap_genes(self, other_genes):
        if self.genome != other_genes.genome:
            raise ValueError(
                "GenomicRegions set-operations only work if both have the same genome. You had %s and %s"
                % (self.genome, other_genes.genome)
            )
        if not "gene_stable_id" in other_genes.df.columns:
            raise ValueError(
                "other_genes must be a Genes (or at least have df[:, 'gene_stable_id']"
            )
        return len(
            set(self.df["gene_stable_id"]).intersection(
                other_genes.df["gene_stable_id"]
            )
        )

    @lazy_method
    def regions_tss(self):
        """Return 'point' regions for the transcription start sites, one per gene"""

        def load():
            res = []
            for dummy_idx, row in self.df.iterrows():
                res.append(
                    {
                        "chr": row["chr"],
                        "start": row["tss"],
                        "stop": row["tss"] + 1,
                        "gene_stable_id": row["gene_stable_id"],
                        "tss_direction": row["strand"],
                    }
                )
            return pd.DataFrame(res)

        return GenomicRegions(
            self.name + " TSS", load, [self.load()], self.genome, on_overlap="merge"
        )

    @lazy_method
    def regions_tes(self):
        """Return 'point' regions for the transcription end sites, one per gene"""

        def load():
            res = []
            for dummy_idx, row in self.df.iterrows():
                res.append(
                    {
                        "chr": row["chr"],
                        "start": row["tes"],
                        "stop": row["tes"] + 1,
                        "gene_stable_id": row["gene_stable_id"],
                    }
                )
            return pd.DataFrame(res)

        return GenomicRegions(
            self.name + " TES", load, [self.load()], self.genome, on_overlap="merge"
        )

    @lazy_method
    def regions_introns(self):
        """Return positions of all intronic regions - possibly overlapping"""

        def load():
            res = {
                "chr": [],
                "start": [],
                "stop": [],
                "transcript_stable_id": [],
                "intron_rank": [],
                "gene_stable_id": [],
            }
            canonical_chromosomes = self.genome.get_chromosome_lengths()
            for (transcript_stable_id, transcript_row) in self.genome.df_transcripts[
                ["gene_stable_id", "chr", "strand"]
            ].iterrows():
                if (
                    not transcript_row["chr"] in canonical_chromosomes
                ):  # pragma: no cover
                    continue
                tr = self.genome.transcripts[transcript_stable_id]
                introns = tr.introns
                cardinality = 0
                for start, stop in introns:
                    res["chr"].append(transcript_row["chr"])
                    res["start"].append(start)
                    res["stop"].append(stop)
                    res["transcript_stable_id"].append(transcript_stable_id)
                    res["gene_stable_id"].append(transcript_row["gene_stable_id"])
                    if transcript_row["strand"] == 1:
                        res["intron_rank"].append(cardinality)
                    else:
                        res["intron_rank"].append(len(introns) - 1 - cardinality)
                    cardinality += 1
            return pd.DataFrame(res)

        return GenomicRegions(
            self.name + " introns",
            load,
            [self.load()],
            self.genome,
            on_overlap="ignore",
        )

    @lazy_method
    def regions_exons_overlapping(self):
        """Return positions of all exonic regions - possibly overlapping"""

        if self.load_strategy.build_deps:
            deps = [
                self.load(),
                ppg.FunctionInvariant(
                    "GenomicRegions_{self.name}_exons_overlapping_actual_load",
                    type(self.genome).df_exons,
                ),
            ]
        else:
            deps = []
        return GenomicRegions(
            self.name + "_exons_overlapping",
            lambda: self.genome.df_exons,
            deps,
            self.genome,
            on_overlap="ignore",
        )

    @lazy_method
    def regions_exons_merged(self):
        """Return positions of all exonic regions - possibly overlapping"""
        if self.load_strategy.build_deps:
            deps = [
                self.load(),
                ppg.FunctionInvariant(
                    "GenomicRegions_{self.name}_exons_merged_actual_load",
                    type(self.genome).df_exons,
                ),
            ]
        else:
            deps = []

        return GenomicRegions(
            self.name + "_exons_merged",
            lambda: self.genome.df_exons,
            deps,
            self.genome,
            on_overlap="merge",
        )

    def write_bed(self, output_filename=None):
        """Store the intervals of the GenomicRegion in a BED file"""
        from mbf_fileformats.bed import BedEntry, write_bed

        output_filename = self.pathify(output_filename, self.name + ".bed")

        def write(output_filename):
            bed_entries = []
            for dummy_idx, row in self.df.iterrows():
                entry = BedEntry(
                    row["chr"], row["start"], row["stop"], name=row["gene_stable_id"]
                )
                bed_entries.append(entry)
            write_bed(output_filename, bed_entries, {}, self.name, include_header=True)

        if self.load_strategy.build_deps:
            deps = [self.load()]
        else:
            deps = []
        return self.load_strategy.generate_file(output_filename, write, deps)

    def write(self, output_filename=None, mangler_function=None):
        if output_filename is None and mangler_function is not None:
            raise ValueError(
                "plain write (without output filename) is always unmangled"
            )
        output_filename = self.pathify(output_filename, self.get_table_filename())
        output_filename.parent.mkdir(exist_ok=True, parents=True)

        def write_mangler_function(df):
            """Establish a decent column order for genes.write"""
            order = [
                "gene_stable_id",
                "name",
                "chr",
                "start",
                "stop",
                "strand",
                "tss",
                "tes",
                "description",
                "biotype",
            ]
            order = [x for x in order if x in df.columns]
            real_order = order + sorted(
                [
                    x
                    for x in df.columns
                    if x not in order
                    and not x.startswith("mysummit")
                    and not x == "parent_row"
                ]
            )
            if "parent_row" in df.columns:
                real_order += ["parent_row"]
            df = df[real_order]
            if mangler_function:
                df = mangler_function(df.copy())
            return df

        return GenomicRegions.write(
            self, output_filename, mangler_function=write_mangler_function
        )
