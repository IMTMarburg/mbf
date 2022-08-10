import pypipegraph as ppg
import pandas as pd
import numpy as np
import math

from .regions import GenomicRegions
from mbf_nested_intervals import merge_df_intervals
from mbf_externals.util import to_string
from ..util import read_pandas


def verify_same_genome(list_of_grs):
    if len(set([x.genome for x in list_of_grs])) > 1:
        raise ValueError("Mixing GenomicRegions from different genomes not supported")


def GenomicRegions_FromGFF(
    name,
    filename,
    genome,
    filter_function=None,
    comment_char=None,
    on_overlap="raise",
    chromosome_mangler=None,
    fix_negative_coordinates=False,
    alternative_class=None,
    summit_annotator=None,
    vid=None,
):
    """Create a GenomicRegions from a gff file.
    You can filter entries with @filter_function(gff_entry_dict) -> Bool,
    remove comment lines starting with a specific character with @comment_char,
    mangle the chromosomes with @chromosome_mangler(str) -> str,
    replace negative coordinates with 0 (@fix_negative_coordinates),
    or provide an alternative constructor to call with @alternative_class
    """

    def load():
        from mbf_fileformats.gff import gffToDict

        entries = gffToDict(filename, comment_char=comment_char)
        data = {
            "chr": [],
            "start": [],
            "stop": [],
            "score": [],
            "strand": [],
            "name": [],
        }
        name_found = False
        for entry in entries:
            if filter_function and not filter_function(entry):
                continue
            if chromosome_mangler:
                chr = chromosome_mangler(entry["seqname"])
            else:
                chr = entry["seqname"]
            data["chr"].append(to_string(chr))
            start = entry["start"]
            if fix_negative_coordinates and start < 0:
                start = 0
            data["start"].append(start)
            data["stop"].append(entry["end"])
            data["score"].append(entry["score"])
            data["strand"].append(entry["strand"])
            name = entry["attributes"].get("Name", [""])[0]
            data["name"].append(name)
            if name:
                name_found = True
        if not name_found:
            del data["name"]
        return pd.DataFrame(data)

    if alternative_class is None:  # pragma: no cover
        alternative_class = GenomicRegions
    if ppg.inside_ppg():
        deps = [
            ppg.FileTimeInvariant(filename),
            ppg.ParameterInvariant(
                name + "_params_GenomicRegions_FromGFF",
                (comment_char, fix_negative_coordinates),
            ),
            ppg.FunctionInvariant(
                name + "_filter_func_GenomicRegions_FromGFF", filter_function
            ),
            ppg.FunctionInvariant(
                name + "_chromosome_manlger_GenomicRegions_FromGFF", chromosome_mangler
            ),
        ]
    else:
        deps = []
    return alternative_class(
        name, load, deps, genome, on_overlap, summit_annotator=summit_annotator, vid=vid
    )


def GenomicRegions_FromWig(
    name,
    filename,
    genome,
    enlarge_5prime=0,
    enlarge_3prime=0,
    on_overlap="raise",
    comment_char=None,
    summit_annotator=None,
    vid=None,
):
    """Create GenomicRegions from a Wiggle file.

    @enlarge_5prime and @enlarge_3prime increase the size of the fragments described in the wig in
    the respective direction (for example if a chip-chip array did not cover every base).
    @comment_char defines which lines to ignore in the wiggle (see {mbf_fileformats.wiggle_to_intervals})

    The resulting GenomicRegions has a column 'Score' that contains the wiggle score"""
    from mbf_fileformats.wiggle import wiggle_to_intervals

    def load():
        df = wiggle_to_intervals(filename, comment_char=comment_char)
        df["chr"] = [to_string(x) for x in df["chr"]]
        df["start"] -= enlarge_5prime
        df["stop"] += enlarge_3prime
        return df

    if ppg.inside_ppg():
        deps = [ppg.FileTimeInvariant(filename)]
    else:
        deps = []

    return GenomicRegions(
        name, load, deps, genome, on_overlap, summit_annotator=summit_annotator, vid=vid
    )


def GenomicRegions_FromBed(
    name,
    filename,
    genome,
    chromosome_mangler=lambda x: x,
    on_overlap="raise",
    filter_invalid_chromosomes=False,
    summit_annotator=None,
    sheet_name=None,
    vid=None,
):
    """Create GenomicRegions from a Bed file.

    The resulting GenomicRegions has a column 'Score' that contains the wiggle score"""
    from mbf_fileformats.bed import read_bed

    def load():
        valid_chromosomes = set(genome.get_chromosome_lengths())
        data = {}
        entries = read_bed(filename)
        data["chr"] = np.array(
            [chromosome_mangler(to_string(e.refseq)) for e in entries], dtype=object
        )
        data["start"] = np.array([e.position for e in entries], dtype=np.int32)
        data["stop"] = np.array(
            [e.position + e.length for e in entries], dtype=np.int32
        )
        data["score"] = np.array([e.score for e in entries], dtype=float)
        data["strand"] = np.array([e.strand for e in entries], dtype=np.int8)
        data["name"] = np.array([to_string(e.name) for e in entries], dtype=object)
        data = pd.DataFrame(data)
        if filter_invalid_chromosomes:  # pragma: no cover
            keep = [x in valid_chromosomes for x in data["chr"]]
            data = data[keep]
        res = data
        if len(res) == 0:
            raise ValueError("Emtpty Bed file - %s" % filename)
        if (np.isnan(res["score"])).all():
            res = res.drop(["score"], axis=1)
        if (len(res["name"]) > 1) and (len(res["name"].unique()) == 1):
            res = res.drop(["name"], axis=1)
        return res

    if ppg.inside_ppg():
        deps = [
            ppg.FileTimeInvariant(filename),
            ppg.FunctionInvariant(name + "_chrmangler", chromosome_mangler),
        ]
    else:
        deps = []

    return GenomicRegions(
        name,
        load,
        deps,
        genome,
        on_overlap=on_overlap,
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
        vid=vid,
    )


def GenomicRegions_FromBigBed(
    name,
    filename,
    genome,
    chromosome_mangler=lambda x: x,
    on_overlap="raise",
    summit_annotator=None,
    sheet_name=None,
    vid=None,
):
    """Create GenomicRegions from a BigBed file.
    @chromosome_mangler translates genome chromosomes into the bigbed's chromosomes!

    """
    from mbf_fileformats.bed import read_bigbed

    def load():
        res = read_bigbed(filename, genome.get_chromosome_lengths(), chromosome_mangler)
        if (res["strand"] == 1).all():
            res = res.drop("strand", axis=1)
        if len(res) == 0:  # pragma: no cover
            raise ValueError(
                "Emtpty BigBed file (or wrong chromosome names)- %s" % filename
            )
        return res

    if ppg.inside_ppg():
        deps = [
            ppg.FileTimeInvariant(filename),
            ppg.FunctionInvariant(name + "_chrmangler", chromosome_mangler),
        ]
    else:
        deps = []

    return GenomicRegions(
        name,
        load,
        deps,
        genome,
        on_overlap=on_overlap,
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
        vid=vid,
    )


def GenomicRegions_BinnedGenome(
    genome, bin_size, limit_to_chromosomes=None, new_name=None, vid=None
):
    """Create GenomicRegions than partition a given chromosome in bins
    of @bin_size bp
    """
    if limit_to_chromosomes and not isinstance(limit_to_chromosomes, list):
        raise ValueError("Limit_to_chromosomes must be a list")

    def load():
        data = {"chr": [], "start": [], "stop": []}
        if limit_to_chromosomes:
            chr_plus_length = [
                (chr, chr_len)
                for (chr, chr_len) in genome.get_chromosome_lengths().items()
                if chr in limit_to_chromosomes
            ]
        else:
            chr_plus_length = genome.get_chromosome_lengths().items()
        if not chr_plus_length:  # pragma: no cover
            raise ValueError("No chromosomes generated")
        for chr, chr_len in chr_plus_length:
            no_of_bins = int(math.ceil(float(chr_len) / bin_size))
            data["chr"].extend([chr] * no_of_bins)
            starts = np.array(range(0, chr_len, bin_size))
            stops = starts + bin_size  # not -1, we are right exclusive!
            data["start"].extend(starts)
            data["stop"].extend(stops)
        res = pd.DataFrame(data)
        return res

    if new_name is None:
        name = genome.name + " binned (%i bp)" % bin_size
    else:
        name = new_name
    return GenomicRegions(name, load, [], genome, vid=vid)


def GenomicRegions_Union(
    name, list_of_grs, summit_annotator=None, sheet_name="Overlaps"
):
    """Combine serveral GRs into one.

    Do not set on_overlap


    """
    verify_same_genome(list_of_grs)

    def load():
        dfs = [x.df[["chr", "start", "stop"]] for x in list_of_grs]
        return pd.concat(dfs, axis=0)

    if ppg.inside_ppg():
        deps = [x.load() for x in list_of_grs]
        deps.append(
            ppg.ParameterInvariant(
                name + "_input_grs", list(sorted([x.name for x in list_of_grs]))
            )
        )
    else:
        deps = []
    vid = ("union", [x.vid for x in list_of_grs])
    return GenomicRegions(
        name,
        load,
        deps,
        list_of_grs[0].genome,
        on_overlap="merge",
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
        vid=vid,
    )


def GenomicRegions_Overlapping(
    new_name, gr_a, gr_b, summit_annotator=None, sheet_name="Overlaps"
):
    return GenomicRegions_Common(new_name, [gr_a, gr_b], summit_annotator, sheet_name)


def GenomicRegions_Common(
    name, list_of_grs, summit_annotator=None, sheet_name="Overlaps"
):
    """Combine serveral GRs into one. Keep only those (union) regions occuring in all."""

    def load():
        union = merge_df_intervals(
            pd.concat([x.df[["chr", "start", "stop"]] for x in list_of_grs])
        ).reset_index(drop=True)
        keep = np.ones((len(union),), dtype=bool)
        for gr in list_of_grs:
            for ii, row in union.iterrows():
                if keep[
                    ii
                ]:  # no point in checking if we already falsified - short circuit...
                    if not gr.has_overlapping(row["chr"], row["start"], row["stop"]):
                        keep[ii] = False
        return union[keep]

    verify_same_genome(list_of_grs)
    if ppg.inside_ppg():
        deps = [x.load() for x in list_of_grs]
        deps.append(
            ppg.ParameterInvariant(
                name + "_input_grs", sorted([x.name for x in list_of_grs])
            )
        )
    else:
        for x in list_of_grs:
            x.load()
        deps = []
    vid = ("common", [x.vid for x in list_of_grs])
    return GenomicRegions(
        name,
        load,
        deps,
        list_of_grs[0].genome,
        on_overlap="raise",
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
        vid=vid,
    )


def GenomicRegions_Invert(new_name, gr, summit_annotator=None, sheet_name="Inverted"):
    """Invert a GenomicRegions. What was covered becomes uncovered, what was uncovered becomes covered.
        [(10, 100), (400, 450)], in a chromosome of size 1000
        becomes
        [(0, 10), (450, 1000)]

        Note that all interval-set based operations (L{union}, L{intersection}, L{difference})
        drop all columns but chr, start, stop (annotators are merged and readded from all sets involved)
        """

    def do_load():
        from mbf_nested_intervals import IntervalSet
        import itertools

        df = gr.df
        joined = []
        for ii, (chr, start, stop) in enumerate(
            zip(df["chr"], df["start"], df["stop"])
        ):
            joined.append(((chr), start, stop, ii))
        joined.sort(key=lambda tup: tup[0])

        out = []
        chr_lengths = gr.genome.get_chromosome_lengths()
        seen = set()
        for chr, sub_group in itertools.groupby(joined, lambda tup: tup[0]):
            args = [x[1:] for x in sub_group]
            iv = IntervalSet.from_tuples_with_id(args)
            new_order = iv.invert(0, chr_lengths[chr]).to_numpy()
            out.append(
                pd.DataFrame({"start": new_order[0], "stop": new_order[1], "chr": chr})
            )
            seen.add(chr)
        for chr in chr_lengths.keys() - seen:
            out.append(
                pd.DataFrame({"start": [0], "stop": [chr_lengths[chr]], "chr": chr})
            )

        return pd.concat(out).reset_index(drop=True)

    if gr.load_strategy.build_deps:
        deps = [
            gr.load(),
            ppg.ParameterInvariant(
                "GenomicRegions_%s_parents" % new_name, (gr.name)
            ),  # so if you swap out the gr, it's detected...
        ]
    else:
        deps = []

    result = GenomicRegions(
        new_name,
        do_load,
        deps,
        gr.genome,
        on_overlap="raise",
        summit_annotator=summit_annotator,
        vid=(["invert"] + gr.vid) if gr.vid is not None else None,
        sheet_name=sheet_name,
    )
    return result


def GenomicRegions_Difference(
    new_name, gr_a, gr_b, summit_annotator=None, sheet_name="difference"
):
    """Create a difference of these intervals wth other_gr's intervalls
        [(10, 100), (400, 450)],
        [(80, 120), (600, 700)]
        becomes
        [(10, 80), (400, 450)]
        (intervals may be split up!)

        Note that all interval-set based operations (L{union}, L{intersection}, L{difference})
        drop all columns but chr, start, stop (annotators are merged and readded from all sets involved)
        """
    verify_same_genome([gr_a, gr_b])

    def do_load():
        new_rows = []
        for idx, row in gr_a.df[["chr", "start", "stop"]].iterrows():
            overlaps = gr_b.get_overlapping(row["chr"], row["start"], row["stop"])
            if not len(overlaps):  # the easy case...
                new_rows.append(row)
            else:
                overlaps = merge_df_intervals(
                    overlaps
                )  # they are now also sorted, so all we need to do is walk them, keep the regions between them (and within the original interval)
                start = row["start"]
                if overlaps.at[0, "start"] <= start and overlaps.at[0, "stop"] > start:
                    start_i = 1
                    start = overlaps.at[0, "stop"]
                else:
                    start_i = 0
                for ii in range(start_i, len(overlaps)):
                    stop = min(overlaps.at[ii, "start"], row["stop"])
                    new_rows.append({"chr": row["chr"], "start": start, "stop": stop})
                    start = overlaps.at[ii, "stop"]
                if start < row["stop"]:
                    new_rows.append(
                        {"chr": row["chr"], "start": start, "stop": row["stop"]}
                    )

                # todo: work the cutting up magic!
                pass
        if new_rows:
            return pd.DataFrame(new_rows)
        else:
            return pd.DataFrame({"chr": [], "start": [], "stop": []})

    if gr_a.load_strategy.build_deps:
        deps = [
            gr_b.load(),
            gr_a.load(),
            ppg.ParameterInvariant(
                "GenomicRegions_%s_parents" % new_name, (gr_a.name, gr_b.name)
            ),  # so if you swap out the gr, it's detected...
        ]
    else:
        gr_b.load()
        deps = []

    result = GenomicRegions(
        new_name,
        do_load,
        deps,
        gr_a.genome,
        on_overlap="merge",
        summit_annotator=summit_annotator,
        vid=["difference"] + list(gr_a.vid) + list(gr_b.vid),
        sheet_name=sheet_name,
    )
    return result


def GenomicRegions_Intersection(
    new_name, gr_a, gr_b, summit_annotator=None, sheet_name="intersection"
):
    """Create an intersection of all intervals...
        [(10, 100), (400, 450)],
        [(80, 120), (600, 700)]
        becomes
        [(80, 100),]

        Note that all interval-set based operations (L{union}, L{intersection}, L{difference})
        drop all columns but chr, start, stop (annotators are merged and readded from all sets involved)
        """
    verify_same_genome([gr_a, gr_b])

    def do_load():
        new_rows = []
        for chr, start, stop in gr_a._iter_intersections(gr_b):
            new_rows.append({"chr": chr, "start": start, "stop": stop})
        if new_rows:
            return pd.DataFrame(new_rows)
        else:
            return pd.DataFrame({"chr": [], "start": [], "stop": []})

    if gr_a.load_strategy.build_deps:
        deps = [
            gr_b.load(),
            gr_a.load(),
            ppg.ParameterInvariant(
                "GenomicRegions_%s_parents" % new_name, (gr_a.name, gr_b.name)
            ),  # so if you swap out the gr, it's detected...
        ]
    else:
        deps = []
        gr_b.load()

    result = GenomicRegions(
        new_name,
        do_load,
        deps,
        gr_a.genome,
        on_overlap="merge",
        summit_annotator=summit_annotator,
        vid=["intersection"] + gr_a.vid + gr_b.vid,
        sheet_name=sheet_name,
    )
    return result


def GenomicRegions_FromPartec(
    name, filename, genome, on_overlap="raise", summit_annotator=None, vid=None
):
    """create GenomicRegions from Partec's output"""
    return GenomicRegions_FromTable(
        name,
        filename,
        genome,
        on_overlap=on_overlap,
        summit_annotator=summit_annotator,
        vid=vid,
        chr_column="Chromosome",
        start_column="Start",
        stop_column="Stop",
        drop_further_columns=False,
        reader=lambda x: pd.read_csv(x, sep="\t"),
    )


def GenomicRegions_FromTable(
    name,
    filename,
    genome,
    on_overlap="raise",
    summit_annotator=None,
    filter_func=None,
    vid=None,
    sheet_name="FromTable",
    drop_further_columns=True,
    chr_column="chr",
    start_column="start",
    stop_column="stop",
    one_based=False,
    reader=read_pandas,
):
    """Read a table file (csv/tsv/xls) with the chr/start/stop columns (renamed?), optionally
    drop all further columns"""

    def load():

        df = reader(filename)
        df["chr"] = df[chr_column].astype(str)
        df["start"] = df[start_column].astype(int)
        if one_based:  # pragma: no cover
            df["start"] -= 1
        df["stop"] = df[stop_column].astype(int)
        if drop_further_columns:  # pragma: no cover
            df = df[["chr", "start", "stop"]]
        if filter_func:  # pragma: no cover
            df = filter_func(df)
        return df

    if ppg.inside_ppg():
        deps = [
            ppg.FileTimeInvariant(filename),
            ppg.FunctionInvariant(name + "_filter_func", filter_func),
        ]
    else:
        deps = []
    return GenomicRegions(
        name,
        load,
        deps,
        genome,
        on_overlap,
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
        vid=vid,
    )


def GenomicRegions_CommonInAtLeastX(
    name, list_of_grs, X, summit_annotator=None, sheet_name="Overlaps"
):
    """Combine serveral GRs into one. Keep only those (union) regions occuring in at least x."""

    def load():
        union = merge_df_intervals(
            pd.concat([x.df[["chr", "start", "stop"]] for x in list_of_grs])
        ).reset_index()
        keep = np.zeros((len(union),), dtype=bool)
        for ii, row in union.iterrows():
            count = 0
            for gr in list_of_grs:
                if gr.has_overlapping(row["chr"], row["start"], row["stop"]):
                    count += 1
            keep[ii] = count >= X
        if not keep.any():  # pragma: no cover
            raise ValueError("Filtered all of them")
        return union.iloc[keep]

    if len(set([x.genome for x in list_of_grs])) > 1:  # pragma: no cover
        raise ValueError("Can only merge GenomicRegions that have the same genome")
    if ppg.inside_ppg():
        deps = [x.load() for x in list_of_grs]
        deps.append(
            ppg.ParameterInvariant(
                name + "_input_grs", sorted([x.name for x in list_of_grs])
            )
        )
    else:
        deps = []
        [x.load() for x in list_of_grs]
    vid = ("common at least %i" % X, [x.vid for x in list_of_grs])
    return GenomicRegions(
        name,
        load,
        deps,
        list_of_grs[0].genome,
        on_overlap="raise",
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
        vid=vid,
    )


def GenomicRegions_Windows(
    genome,
    name,
    window_size,
    window_spacing,
    subset_of_chromosomes=None,
    sheet_name="Windowed",
):
    """Create a GenomicRegions that has a window of size @window_size (0 for next to each other), windows are spaced @window_spacing
    across all or a @subset_of_chromosomes"""
    if isinstance(genome, str):  # pragma: no cover
        raise TypeError("Wrong arguments, genome must be first")

    def load():
        chrs_to_include = list(
            subset_of_chromosomes
            if subset_of_chromosomes
            else genome.get_chromosome_lengths().keys()
        )
        chrs = []
        starts = []
        stops = []
        chr_lengths = genome.get_chromosome_lengths()
        for c in chrs_to_include:
            ll = chr_lengths[c]
            for ii in range(0, ll, window_size + window_spacing):
                chrs.append(c)
                starts.append(ii)
                stops.append(ii + window_size)
        return pd.DataFrame({"chr": chrs, "start": starts, "stop": stops})

    return GenomicRegions(
        name, load, [], genome, on_overlap="raise", sheet_name=sheet_name
    )


def GenomicRegions_FilterRemoveOverlapping(
    new_name, gr_a, other_grs, summit_annotator=None, sheet_name="Overlaps"
):
    """Filter all from this GenomicRegions that have an overlapping region in other_gr
    Note that filtering does not change the coordinates, it only filters,
    non annotator additional rows are kept, annotators are recalculated.
    """
    if isinstance(other_grs, GenomicRegions):
        other_grs = [other_grs]
    verify_same_genome([gr_a] + other_grs)

    def filter_func(df):
        keep = np.zeros((len(df)), dtype=bool)
        for ii, tup in enumerate(df[["chr", "start", "stop"]].itertuples()):
            for other_gr in other_grs:
                keep[ii] = keep[ii] | other_gr.has_overlapping(
                    tup.chr, tup.start, tup.stop
                )
        return ~keep

    if not summit_annotator:  # pragma: no cover
        summit_annotator = gr_a.summit_annotator

    if gr_a.load_strategy.build_deps:
        deps = [
            [x.load() for x in other_grs],
            ppg.ParameterInvariant(
                "GenomicRegions_%s_parents" % new_name,
                (gr_a.name, [x.name for x in other_grs]),
            ),  # so if you swap out the gr, it's detected...
        ]
    else:
        for other_gr in other_grs:
            other_gr.load()
        deps = []

    return gr_a.filter(
        new_name,
        df_filter_function=filter_func,
        dependencies=deps,
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
    )


def GenomicRegions_FilterToOverlapping(
    new_name, gr_a, other_grs, summit_annotator=None, sheet_name="Overlaps"
):
    """Filter to just those that overlap one in *all* other_grs.
    Note that filtering does not change the coordinates, it only filters,
    non annotator additional rows are kept, annotators are recalculated.
    """
    if isinstance(other_grs, GenomicRegions):
        other_grs = [other_grs]
    verify_same_genome([gr_a] + other_grs)

    def filter_func(df):
        keep = np.ones((len(df)), dtype=bool)
        for ii, row in enumerate(df[["chr", "start", "stop"]].itertuples()):
            for gr in other_grs:
                keep[ii] &= gr.has_overlapping(row.chr, row.start, row.stop)
        return keep

    if gr_a.load_strategy.build_deps:
        deps = [other_gr.load() for other_gr in other_grs] + [
            ppg.ParameterInvariant(
                "GenomicRegions_%s_parents" % new_name,
                (gr_a.name, [other_gr.name for other_gr in other_grs]),
            )  # so if you swap out the gr, it's detected...
        ]
    else:
        for other_gr in other_grs:
            other_gr.load()
        deps = []

    return gr_a.filter(
        new_name,
        df_filter_function=filter_func,
        dependencies=deps,
        summit_annotator=summit_annotator,
        sheet_name=sheet_name,
    )
