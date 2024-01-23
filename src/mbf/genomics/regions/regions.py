import pypipegraph as ppg
import numpy as np
import pandas as pd
from pathlib import Path
import itertools
from collections.abc import Iterator


from mbf.genomics.delayeddataframe import DelayedDataFrame
from mbf.genomes import GenomeBase
from mbf.externals.util import lazy_property
from mbf_nested_intervals import (
    merge_df_intervals,
    merge_df_intervals_with_callback,
    IntervalSet,
    _df_to_tup_no_strand,
)


def merge_identical_but_raise_on_further_overlap(df):
    def iv_func(iv):
        iv = iv.remove_duplicates()
        if iv.any_overlapping():
            raise ValueError(
                "Overlapping after merging identical present. Use  a different on_overlap or fix our data"
            )
        return iv

    res = merge_df_intervals(df, iv_func)
    return res


region_registry = {}


class GenomicRegions(DelayedDataFrame):

    """A genomic regions object encapsulates intervals across the chromosomes
    of a genome. They model everything from restriction sites over binding sites
    to gene locations.

    Internally a genomic region consists of a dataframe {chr, start, stop} (.df),
    and GR-annotators can be used to add further columns.
    Please note that the intervals are 'standard python' ie. left inclusive, right exclusive.
    An interval 10:20, is 10 bp long, and covers 10, 11, ..., 18, 19.
    That also means a 'point' interval is actually x, x+1, not x,x (that's an empty interval)

    There is a variety of set operations you can perform on them that are done efficiently
    via interval trees.

    GRs do not require disjoint regions. More specifically, they can be configured
    to do a variety of things if they're passed overlapping regions - ignore them,
    raise an exception, merge them...

    You can query them for all of their internal regions overlapping an interval

    Annotators are inherited by descendands, their values are recalculated for
    the new intervals.

    """

    # basic loading
    def __init__(
        self,
        name,
        loading_function,
        dependencies,
        genome,
        on_overlap="raise",
        result_dir=None,
        summit_annotator=None,
        sheet_name=None,
        vid=None,
    ):
        """Create a lazy loaded GenomicRegion, available once it's load() job has been completed.
        @dependencies allow you to inject dependencies for the load() job.
        @on_overlap may be one of (raise, ignore, merge, drop, see below), which decides how overlapping regions are handled.

        @loading_function is a parameterless function that returns dataframe of {chr, start, stop}.
        The dataframe does not need to sort the data in any way. chr must be a string column, start and stop integers
        (this is checked).

        The internal df is always sorted by chr, start...


        @on_overlap, allowed values:
            -raise: - raise a ValueError if any overlapping regions are detected (default)
            -merge: - combine overlapping regions into larger ones
            -merge_identical: combine identical regions, but raise if there are non identical overlapping ones
            -(merge, merge_function) - combine regions - call merge_function(sub_df) which must return a row (dict) to keep and is only called when there are multiple rows to pick! Its row['start'], row['stop'] are ignored though
            -drop: - drop all regions that are overlapping one another
            -ignore: - ignore overlapping regions (this implies a nested list search for overlap queries
            and will lead to certain functions that assume non-overlapping regions (such as covered_bases)
            to raise ValueErrors)

        """
        if (
            isinstance(dependencies, str)
            or not hasattr(dependencies, "__iter__")
            or (
                hasattr(dependencies, "__iter__") and isinstance(dependencies, Iterator)
            )
        ):
            raise ValueError(
                "dependencies must be iterable (use [] for no external dependencies"
            )
        if not isinstance(dependencies, list):
            dependencies = list(dependencies)

        allowed_overlap_modes = ("raise", "merge", "ignore", "drop", "merge_identical")
        if not on_overlap in allowed_overlap_modes and not (
            isinstance(on_overlap, tuple)
            and on_overlap[0] == "merge"
            and hasattr(on_overlap[1], "__call__")
        ):
            raise ValueError(
                "Invalid on_overlap mode %s. Allowed: %s, or a tuple of ('merge', function)"
                % (on_overlap, allowed_overlap_modes)
            )

        self.name = name
        self.gr_loading_function = loading_function
        self.genome = genome
        self.on_overlap = on_overlap
        self._default_mangler = True
        if self.on_overlap == "ignore":
            self.need_to_handle_overlapping_regions = True
        elif (
            self.on_overlap == "raise"
            or self.on_overlap == "merge"
            or self.on_overlap == "merge_identical"
            or self.on_overlap == "drop"
            or (
                isinstance(self.on_overlap, tuple)
                and self.on_overlap[0] == "merge"
                and hasattr(self.on_overlap[1], "__call__")
            )
        ):
            self.need_to_handle_overlapping_regions = False
        else:  # pragma: no cover - defensive
            raise ValueError(
                "Don't know how to decide on has_overlapping from %s"
                % (self.on_overlap,)
            )

        if result_dir:
            result_dir = Path(result_dir)
        else:
            if sheet_name:
                result_dir = Path("results") / "GenomicRegions" / sheet_name / name
            else:
                result_dir = Path("results") / "GenomicRegions" / name
        self.sheet_name = sheet_name

        super().__init__(name, self._load, dependencies, result_dir)
        if self.load_strategy.build_deps:
            dependencies = []
            if isinstance(on_overlap, tuple):
                dependencies.append(
                    ppg.ParameterInvariant("grload_params" + name, ("function",))
                )
                dependencies.append(
                    ppg.FunctionInvariant("grload_overlap_func" + name, on_overlap[1])
                )
            else:
                dependencies.append(
                    ppg.ParameterInvariant("grload_params" + name, (on_overlap,))
                )
            dependencies.extend(genome.download_genome())
            dependencies.append(
                ppg.FunctionInvariant("grload_func" + name, self.gr_loading_function)
            )
            self.load_strategy.load().depends_on(dependencies)

        if not hasattr(self, "column_properties"):
            self.column_properties = {
                "chr": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "On which chromosome (or contig) the region is loacted",
                },
                "start": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "Left most position of this region",
                },
                "stop": {
                    "user_visible": False,
                    "priority": -997,
                    "description": "Right most position of this region",
                },
            }
        # self.genome = genome

        self.random_count = 0
        if summit_annotator:
            self.summit_annotator = summit_annotator
        elif summit_annotator is False:
            self.summit_annotator = None
        else:
            from .annotators import SummitMiddle

            self.summit_annotator = SummitMiddle()
        if self.summit_annotator is not None:
            self.add_annotator(self.summit_annotator)
        self.register()
        if vid is None:
            vid = []
        self.vid = vid

    def get_default_columns(self):
        return ("chr", "start", "stop")

    def _load(self):

        df = self.gr_loading_function()
        if not isinstance(df, pd.DataFrame):
            raise ValueError(
                "GenomicRegion(%s).loading_function must return a pandas.DataFrame, was: %s\n%s"
                % (self.name, type(df), self.gr_loading_function)
            )
        for col in self.get_default_columns():
            if not col in df.columns:
                raise ValueError(
                    "%s not in dataframe returned by GenomicRegion.loading_function"
                    % col
                )
        allowed_chromosomes = set(self.genome.get_chromosome_lengths().keys())
        if len(df):
            # if isinstance(df.iloc[0]["chr"], six.string_types):
            # df["chr"] = [str(x) for x in df["chr"].values]
            if not (isinstance(df.iloc[0]["chr"], str)):
                raise ValueError(
                    "Chromosomes must be a string, was: %s, first instance %s of type %s"
                    % (
                        df["chr"].dtype,
                        repr(df.iloc[0]["chr"]),
                        type(df.iloc[0]["chr"]),
                    )
                )
            chrs = set([str(x) for x in df["chr"].unique()])
            if chrs.difference(allowed_chromosomes):
                raise ValueError(
                    "Invalid chromosomes found: %s, expected one of: %s"
                    % (
                        chrs.difference(allowed_chromosomes),
                        sorted(allowed_chromosomes),
                    )
                )
            if not np.issubdtype(df["start"].dtype, np.integer):
                raise ValueError(
                    "start needs to be an integer, was: %s" % df["start"].dtype
                )
            if not np.issubdtype(df["stop"].dtype, np.integer):
                raise ValueError(
                    "stop needs to be an integer, was: %s" % df["stop"].dtype
                )
            if (df["start"] < 0).any():
                print(df[df["start"] < 0])
                raise ValueError("All starts need to be positive!")
            if (df["start"] > df["stop"]).any():
                error = df[(df["start"] > df["stop"])][["chr", "start", "stop"]]
                error = error[:100]  # not too many...
                error = str(error)
                raise ValueError(
                    "%s.loading_function returned a negative interval:\n%s"
                    % (self, error)
                )
            df = self.handle_overlap(df)
        else:
            df = df.assign(
                start=np.array([], dtype=np.int32),
                stop=np.array([], dtype=np.int32),
                chr=np.array([], dtype=object),
            )
        # enforce column order
        cols = ["chr", "start", "stop"]
        for x in df.columns:
            if not x in cols:
                cols.append(x)
        df = df[cols]

        return df

    def handle_overlap(self, df):
        """depending on L{GenomicRegion.on_overlap}, check
        for overlapping regions and handle accordingly"""
        if self.on_overlap == "raise":
            df = df.sort_values(["chr", "start"], ascending=[True, True])
            last_chr = None
            last_stop = 0
            last_row = None
            # todo: vectorize this... even with a stupid group by we should be able to 
            # do this much faster in the ok path...
            for idx, row in df.iterrows():
                if row["chr"] != last_chr:
                    last_chr = row["chr"]
                    last_stop = 0
                if row["start"] < last_stop:
                    raise ValueError(
                        "%s: Overlapping intervals: %s %i..%i vs %i..%i"
                        % (
                            self.name,
                            row["chr"],
                            row["start"],
                            row["stop"],
                            last_row["start"],
                            last_row["stop"],
                        )
                    )
                last_stop = row["stop"]
                last_row = row
            return df
        elif self.on_overlap == "merge":
            return merge_df_intervals(df)
        elif self.on_overlap == "merge_identical":
            return merge_identical_but_raise_on_further_overlap(df)
        elif (
            isinstance(self.on_overlap, tuple)
            and self.on_overlap[0] == "merge"
            and hasattr(self.on_overlap[1], "__call__")
        ):
            return merge_df_intervals_with_callback(df, self.on_overlap[1])

        elif self.on_overlap == "ignore":
            df = df.sort_values(["chr", "start"], ascending=[True, True])
            ov = []
            tups = _df_to_tup_no_strand(df)
            for chr, sub_group in itertools.groupby(tups, lambda tup: tup[0]):
                args = [x[1:] for x in sub_group]
                for r in args:
                    if r[0] > r[1]:
                        print('error tuple!', r)

                try:
                    iv = IntervalSet.from_tuples_with_id(args)
                except:
                    raise
                ov.extend(iv.overlap_status())
            df = df.assign(is_overlapping=ov)
            return df

        elif self.on_overlap == "drop":
            return merge_df_intervals(df, lambda iv: iv.merge_drop())
        else:  # pragma: no branch - defensive
            raise NotImplementedError(  # pragma: no cover
                "This branch should not happen - unhandled on_overlap"
            )

    def do_build_intervals(self):
        """ "Build the interval trees right now, ignoring all dependencies"""
        if not hasattr(self, "_interval_sets"):
            self._interval_sets = {}
            for chr, tups in itertools.groupby(
                self.df[["chr", "start", "stop"]].itertuples(), key=lambda tup: tup.chr
            ):
                tups = [(tup.start, tup.stop, tup[0]) for tup in tups]
                self._interval_sets[chr] = IntervalSet.from_tuples_with_id(tups)

    def has_overlapping(self, chr, start, stop):
        """is there an interval overlapping the region passed"""
        self.do_build_intervals()
        if not chr in self._interval_sets:
            return False
        return self._interval_sets[chr].has_overlap(start, stop)

    def get_overlapping(self, chr, start, stop):
        """Retrieve the rows for the region passed in.
        Returns an empty dataframe if there is no overlap

        Please note that the interval is end excluding - ie start == stop
        means an empty interval and nothing ever overlapping!
        """
        # print 'testing overlap for', chr, start, stop
        self.do_build_intervals()
        if not chr in self._interval_sets:
            return self.df[0:0]
        ids = sorted(self._interval_sets[chr].get_overlap(start, stop).to_ids())
        return self.df.loc[ids]

    def get_closest_by_start(self, chr, point):
        """Find the interval that has the closest to the passed point.
        Returns a df with that interval, or an empty df!
        """
        # first we check whether we're on top of a region...
        # overlapping = self.get_overlapping(chr, point, point+1)
        # if len(overlapping):
        # return overlapping # could be more than one?!
        self.do_build_intervals()
        if not chr in self._interval_sets:
            return self.df[0:0]
        res = self._interval_sets[chr].find_closest_start(point)
        if res is None:  # pragma: no cover
            return self.df[0:0]
        start, stop, ids = res
        return self.df.loc[ids]

    # various statistics
    def get_no_of_entries(self):
        """How many intervals are there"""
        return len(self.df)

    @lazy_property
    def covered_bases(self):
        """How many base pairs are covered by these intervals"""
        if self.on_overlap == "ignore":
            raise ValueError(
                "covered_bases is currently only implemented for not overlapping GenomicRegions (on_overlap != 'ignore')"
            )
        return (self.df["stop"] - self.df["start"]).sum()

    @lazy_property
    def mean_size(self):
        """Get the mean size of the intervals defined in this GR"""
        return np.mean(self.df["stop"] - self.df["start"])

    def register(self):
        region_registry[self.name] = self

    # magic
    def __hash__(self):
        return hash("GR" + self.name)

    def __str__(self):
        return "GenomicRegion(%s)" % self.name

    def __repr__(self):
        return "GenomicRegion(%s)" % self.name

    # output functions

    def write_bed(self, output_filename=None, region_name=None, include_header=False):
        """Store the intervals of the GenomicRegion in a BED file
        @region_name: Insert the column of the GenomicRegions-Object e.g. 'repeat_name'"""
        from mbf.fileformats.bed import BedEntry, write_bed

        output_filename = self.pathify(output_filename, self.name + ".bed").relative_to(
            Path(".").absolute()
        )

        def write(output_filename):
            bed_entries = []
            for ii, row in self.df.iterrows():
                if region_name is None:

                    def get_name(row):
                        return None

                else:

                    def get_name(row):
                        if region_name in row:
                            return row[region_name]
                        else:
                            raise KeyError(
                                "key: %s not in genomic regions object. These objects are available: %s"
                                % (region_name, self.df.columns)
                            )

                entry = BedEntry(
                    row["chr"],
                    row["start"],
                    row["stop"],
                    name=get_name(row),
                    strand=row["strand"] if "strand" in row else None,
                    score=row["score"] if "score" in row else None,
                )
                bed_entries.append(entry)
            write_bed(
                output_filename,
                bed_entries,
                {},
                self.name,
                include_header=include_header,
            )

        if self.load_strategy.build_deps:
            deps = [
                self.load(),
                ppg.ParameterInvariant(output_filename, (include_header,)),
            ]
        else:
            deps = []
        return self.load_strategy.generate_file(output_filename, write, deps)

    def write_bigbed(
        self, output_filename=None, name_column=None
    ):  # pragma: no cover - till we have the test
        """Store the intervals of the GenomicRegion in a big bed file"""
        from mbf.fileformats.bed import BedEntry, write_bigbed

        output_filename = self.pathify(
            output_filename, self.name + ".bigbed"
        ).relative_to(Path(".").absolute())

        def write(output_filename):

            bed_entries = []
            for idx, row in self.df.iterrows():
                entry = BedEntry(
                    row["chr"],
                    row["start"],
                    row["stop"],
                    name=row[name_column] if name_column else None,
                    strand=row["strand"] if "strand" in row else 0,
                )
                bed_entries.append(entry)
            if len(self.df):
                write_bigbed(
                    bed_entries, output_filename, self.genome.get_chromosome_lengths()
                )
            else:
                with open(output_filename, "wb"):
                    pass

        if self.load_strategy.build_deps:
            deps = [self.load()]
        else:
            deps = []
        return self.load_strategy.generate_file(
            output_filename, write, deps, empty_ok=True
        )

    # filtering
    def _new_for_filtering(self, new_name, load_func, deps, **kwargs):
        """When filtering, a new object of this class is created.
        To pass it the right options from the parent, overwrite this
        """
        kwargs["sheet_name"] = kwargs.get("sheet_name", "Filtered")
        for k in "on_overlap", "result_dir", "summit_annotator", "vid":
            if not k in kwargs:
                kwargs[k] = getattr(self, k)
        return GenomicRegions(new_name, load_func, deps, genome=self.genome, **kwargs)

    def _iter_intersections(self, other_gr):
        """Iterate over (chr, start, stop) tuples of the intersections between this GenomicRegions
        and other_gr.

        Refactored from intersection and overlap_basepairs"""
        if self.genome != other_gr.genome:
            raise ValueError(
                "GenomicRegions set-operations only work if both have the same genome. You had %s and %s"
                % (self.genome, other_gr.genome)
            )
        if self.need_to_handle_overlapping_regions:
            raise ValueError(
                "_iter_intersections currently only works on non_overlapping intervalsets - %s was %s"
                % (self.name, self.on_overlap)
            )
        if other_gr.need_to_handle_overlapping_regions:
            other_df = merge_df_intervals(other_gr.df)
        else:
            other_df = other_gr.df
        for chr in self.df["chr"].unique():
            df_here = self.df[self.df["chr"] == chr]
            starts_here = df_here["start"]
            stops_here = df_here["stop"]
            print(other_df)
            print(other_df.dtypes)
            df_there = other_df[other_df["chr"] == chr]
            starts_there = df_there["start"]
            stops_there = df_there["stop"]
            ii = 0
            jj = 0
            while ii < len(df_here) and jj < len(df_there):
                if starts_there.iloc[jj] <= starts_here.iloc[ii] < stops_there.iloc[jj]:
                    yield (
                        chr,
                        max(starts_here.iloc[ii], starts_there.iloc[jj]),
                        min(stops_here.iloc[ii], stops_there.iloc[jj]),
                    )
                    ii += 1
                elif (
                    starts_here.iloc[ii] <= starts_there.iloc[jj] <= stops_here.iloc[ii]
                ):
                    yield (
                        chr,
                        max(starts_here.iloc[ii], starts_there.iloc[jj]),
                        min(stops_here.iloc[ii], stops_there.iloc[jj]),
                    )
                    jj += 1
                elif starts_here.iloc[ii] < starts_there.iloc[jj]:
                    ii += 1
                else:
                    jj += 1
        return

    def copy_annotators(self, *other_grs):
        """Copy annotators from other GRs"""
        for p in other_grs:
            for anno in p.annotators.values():
                self += anno

    # overlap calculations
    def overlap_basepair(
        self, other_gr
    ):  # todo: given to GRs with on_overlap != ignore, we could do a sorted search like we did for the GIS, that should be faster...
        """calculate the overlap between to GenomicRegions, on a base pair level - ie. the size of the intersection set
        (this converts other_gr into a disjoint set, so each base is counted only once. self must be disjoint)"""
        if self.genome != other_gr.genome:
            raise ValueError(
                "GenomicRegions set-operations only work if both have the same genome. You had %s and %s"
                % (self.genome, other_gr.genome)
            )
        overlap = 0
        for chr, start, stop in self._iter_intersections(other_gr):
            overlap += stop - start
        return overlap

    def overlap_percentage(self, other_gr):
        """calculate the percentage of overlapping base pairs
        Basically overlap in bp / min(bp self, bp other_gr)
        if either is empty, percentage = 0
        """
        if self.genome != other_gr.genome:
            raise ValueError(
                "GenomicRegions set-operations only work if both have the same genome. You had %s and %s"
                % (self.genome, other_gr.genome)
            )
        bp_overlap = self.overlap_basepair(other_gr)
        divisor = min(self.covered_bases, other_gr.covered_bases)
        if divisor == 0:
            return 0
        else:
            return float(bp_overlap) / divisor

    def overlap_count(self, other_gr):
        """Count the size of the union of all intersecting intervals. Basically, len(overlapping(..., other_gr))

        See L{overlapping}.

        Requires self.load() to have been done.

        This is (much) faster than building overlap()  and calling len,.

        The invariant len(a) + len(b) == len(a.filter_remove_overlapping(b))
            + len(b.filter_remove_overlapping(a)) + a.overlap_count(b) is not true,
             because the overlapping count only counts each region once,
             even if it was compromised out of multiple regions in the input. Say you have
        ---  ----  -----
          ----  ----
        that count's as 3 (or two) regions in len(a) (len(b)),
         but just as one in overlapping/overlap_count (and the filtered ones are of len(0)). Still, I argue that this is the right one to use in venn diagrams.

        """
        if self.genome != other_gr.genome:
            raise ValueError(
                "GenomicRegions set-operations only work if both have the same genome. You had %s and %s"
                % (self.genome, other_gr.genome)
            )
        overlap = 0
        for chrom in self.genome.get_chromosome_lengths():
            a = self.df[self.df.chr == chrom][["start", "stop"]].assign(color=1)
            b = other_gr.df[other_gr.df.chr == chrom][["start", "stop"]].assign(color=2)
            comb = pd.concat([a, b]).sort_values("start")
            last_stop = 0
            run_start = -1
            run_colors = 0
            for dummy_idx, row in comb.iterrows():
                start = row["start"]
                stop = row["stop"]
                color = row["color"]
                if start < last_stop:  # a continuation
                    run_colors |= color
                    last_stop = max(last_stop, stop)
                else:  # end of this run
                    if run_start > 0:
                        if run_colors == 3:
                            overlap += 1
                    run_start = start
                    last_stop = stop
                    run_colors = color
            if (
                len(comb) and start != run_start and run_colors == 3
            ):  # capture the last one that will not have triggered an 'end-of-run'
                overlap += 1
        return overlap

    def intersection_count(self, other_gr):
        """calculate the number of intersecting regions between two GenomicRegions
        Note that
        len(here) - intersection_count(here, other_gr) != len(here.filter_remove_overlapping(...,other_gr)
        since they handle this situation differently:
        here : ---    ----
        there:   ------
        (which here & in intersection gives you two intervals, while in overlapping (and overlap_count) it
        gives just one)

        """
        if self.genome != other_gr.genome:
            raise ValueError(
                "GenomicRegions set-operations only work if both have the same genome. You had %s and %s"
                % (self.genome, other_gr.genome)
            )
        if not self.load_strategy.build_deps:
            other_gr.load()
        overlap = 0
        for chr, start, stop in self._iter_intersections(other_gr):
            overlap += 1
        return overlap

    def convert(
        self,
        new_name,
        conversion_function,
        new_genome=None,
        dependencies=None,
        on_overlap="raise",
        sheet_name="Converted",
        anno_columns_to_keep=[],
        summit_annotator=None,
    ):
        """Convert the intervals into some other coordinate system,
        possibly changing the genome as well.
        Does not conserve annotators, non canonical rows depend on the conversion fucntion

        @conversion_function is passed the dataframe, and must return one containing
        at least (chr, start, stop).
        if conversion_func.take_genomicregion is True, it get's passed the GenomicRegion instead. 
          (it must still return a dataframe)
        @conversion_function may be a tuple (function, [annotators]),
            in which case the function is treated as conversion_function
            and the annotators are added as dependencies
        @conversion_function may be a tuple (function, [annotators], parameters),
            in which case the function is treated as conversion_function and the
            annotators are added as dependencies,
            and an additional parameter_dependency is added
        @dependencies must be a list of jobs
        """
        if not isinstance(new_name, str):
            raise ValueError("Name must be a string")
        if isinstance(conversion_function, tuple):
            if len(conversion_function) > 2:
                convert_parameters = conversion_function[2]
            else:
                convert_parameters = None
            annotators_required = conversion_function[1]
            conversion_function = conversion_function[0]
        else:
            annotators_required = []
            convert_parameters = None

        def do_load():
            if getattr(conversion_function, 'take_genomicregion', False):
                df = conversion_function(self)
            else:
                df = conversion_function(self.df)
            if not isinstance(df, pd.DataFrame):
                raise ValueError(
                    "GenomicRegions.convert conversion_function must return a pandas.DataFrame."
                )
            for col in df.columns[:]:
                if not (
                    col in self.non_annotator_columns or col in anno_columns_to_keep
                ):
                    df = df.drop(col, axis=1)
            return df

        if new_genome is None:
            new_genome = self.genome
        else:
            if not isinstance(new_genome, GenomeBase):
                raise ValueError(
                    "new_genome %s was not a genome. Did you mean to pass in dependencies, that's the next parameter?"
                    % new_genome
                )
        if self.load_strategy.build_deps:
            deps = [self.load()]
            if dependencies:
                deps.extend(dependencies)
            for anno in annotators_required:
                deps.append(self.add_annotator(anno))
            deps.append(
                ppg.ParameterInvariant(
                    new_name + "_conversion_paramaters", convert_parameters
                )
            )
            deps.append(
                ppg.FunctionInvariant(
                    new_name + "_conversion_function", conversion_function
                )
            )
        else:
            deps = []
        if hasattr(conversion_function, "dependencies"):
            deps.extend(conversion_function.dependencies)
        if summit_annotator is None:  # pragma: no cover
            summit_annotator = self.summit_annotator
        return GenomicRegions(
            new_name,
            do_load,
            deps,
            new_genome,
            on_overlap=on_overlap,
            summit_annotator=summit_annotator,
            sheet_name=sheet_name,
            vid=self.vid,
        )
