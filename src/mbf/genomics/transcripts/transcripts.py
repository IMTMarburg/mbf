from ..regions import GenomicRegions
import numpy as np
import pypipegraph as ppg
import pandas as pd
from pathlib import Path
from mbf_externals.util import lazy_method

_transcripts_per_genome_singletons = {}


class Transcripts(GenomicRegions):
    def __new__(cls, genome, alternative_load_func=None, *args, **kwargs):
        """Make sure that Genes for a full genome (ie. before filtering) are singletonic. Ie.
        you can always safely call Genes(my_genome), without worrying about duplicate objects"""
        if alternative_load_func is None:
            if ppg.util.global_pipegraph:
                if not hasattr(
                    ppg.util.global_pipegraph, "_transcripts_per_genome_singletons"
                ):
                    ppg.util.global_pipegraph._transcripts_per_genome_singletons = {}
                singleton_dict = (
                    ppg.util.global_pipegraph._transcripts_per_genome_singletons
                )
            else:
                singleton_dict = _transcripts_per_genome_singletons
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
                and alternative_load_func != self.transcripts_load_func
            ):  # pragma: no cover -
                # this can only happen iff somebody starts to singletonize the Transcript
                # with loading_functions - otherwise the duplicate object
                # checker will kick in first
                raise ValueError(
                    "Trying to define Transcript(%s) twice with different loading funcs"
                    % self.name
                )
            pass
        else:
            if name is None:
                if alternative_load_func:
                    raise ValueError(
                        "If you pass in an alternative_load_func you also need to specify a name"
                    )
                name = "Transcript%s" % genome.name
                self.top_level = True
            else:
                self.top_level = False
            if alternative_load_func is None:
                load_func = lambda: genome.df_transcripts.reset_index()  # noqa: E731
            else:
                load_func = alternative_load_func

            self.transcripts_load_func = load_func
            if result_dir:
                pass
            elif sheet_name:
                result_dir = (
                    result_dir or Path("results") / "Transcripts" / sheet_name / name
                )
            else:
                result_dir = result_dir or Path("results") / "Transcripts" / name
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
                "transcript_stable_id": {
                    "user_visible": True,
                    "priority": -1000,
                    "index": True,
                    "description": "Unique identification name for this transcript",
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
        return "Transcripts(%s)" % self.name

    def __repr__(self):
        return "Transcripts(%s)" % self.name

    def register(self):
        pass

    def _load(self):
        """Load func"""
        if hasattr(self, "df"):  # pragma: no cover
            return
        df = self.transcripts_load_func()

        if not isinstance(df, pd.DataFrame):
            raise ValueError(
                "GenomicRegion.loading_function must return a DataFrame, was: %s"
                % type(df)
            )
        if not 'tss' in df.columns and not 'tes' in df.columns and 'start' in df.columns and 'stop' in df.columns and 'strand' in df.columns:
            tss = pd.Series(-1, index=df.index)
            tes = pd.Series(-1, index=df.index)
            forward = df.strand == 1
            tss[forward] = df['start'][forward]
            tss[~forward] = df['stop'][~forward]
            tes[forward] = df['stop'][forward]
            tes[~forward] = df['start'][~forward]
            df = df.assign(tss = tss, tes = tes)

        for col in self.get_default_columns():
            if not col in df.columns:
                func_filename = self.transcripts_load_func.__code__.co_filename
                func_line_no = self.transcripts_load_func.__code__.co_firstlineno
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
                            self.transcripts_load_func,
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
            else:
                assert (df['start'] == np.min([df["tss"], df["tes"]], 0)).all()
            if not "stop" in df.columns:
                df = df.assign(stop=np.max([df["tss"], df["tes"]], 0))
        else:
            df = df.assign(
                start=np.array([], dtype=np.int32), stop=np.array([], dtype=np.int32)
            )

        if (df["start"] > df["stop"]).any():
            raise ValueError(
                "Transcripts.loading_function returned a negative interval:\n %s"
                % df[df["start"] > df["stop"]].head()
            )
        self.df = df.sort_values(["chr", "start"], ascending=[True, True]).reset_index(
            drop=True
        )  # since we don't call handle_overlap
        # enforce column order
        cols = [
            "gene_stable_id",
            "transcript_stable_id",
            "name",
            "chr",
            "start",
            "stop",
            "strand",
            "tss",
            "tes",
        ]
        for x in df.columns:
            if not x in cols:
                cols.append(x)
        df = df[cols]
        return df

    def get_default_columns(self):
        return ("chr", "tss", "tes", "gene_stable_id", "transcript_stable_id", "name")

    def _new_for_filtering(self, new_name, load_func, dependencies, **kwargs):
        """When filtering, a new object of this class is created.
        To pass it the right options from the parent, overwrite this
        """
        return Transcripts(self.genome, load_func, new_name, dependencies, **kwargs)

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
                        "transcript_stable_id": row["transcript_stable_id"],
                        "tss_direction": row["strand"],
                    }
                )
            return pd.DataFrame(res)

        return GenomicRegions(
            self.name + " TSS", load, [self.load()], self.genome, on_overlap="ignore"
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
                        "transcript_stable_id": row["transcript_stable_id"],
                    }
                )
            return pd.DataFrame(res)

        return GenomicRegions(
            self.name + " TES", load, [self.load()], self.genome, on_overlap="ignore"
        )
