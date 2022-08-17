"""Utility functions to pass to GenomicRegions.convert(..., convert_func)

"""

# TODO: Liftover utility -> mbf.externals
# TODO: chains not included?

import tempfile
import pandas as pd
import pypipegraph as ppg
import numpy as np
import subprocess
from mbf.externals.util import to_string, to_bytes
from pathlib import Path

file_path = Path(__file__).parent


def grow(basepairs):
    """A function for GenomicRegions.convert that enlarges the regions
    bei doing start = start - basepairs, stop = stop + basepairs"""

    def do_grow(df):
        starts = df["start"] - basepairs
        starts[starts < 0] = 0
        stops = df["stop"] + basepairs
        new_df = df.copy()
        new_df = new_df.assign(start=starts, stop=stops)
        return new_df

    return do_grow, [], basepairs


def promotorize(upstream=1250, downstream=0):
    """Genes.convert - returns [-upstream...tss+downstream] regions"""

    def do_promotorize(df):
        res = {"chr": df["chr"]}
        res["start"] = np.zeros((len(df),), dtype=np.int32)
        res["stop"] = np.zeros((len(df),), dtype=np.int32)
        forward = df["strand"] == 1
        res["start"][:] = df["tss"]  # Assign within array.
        res["stop"][:] = df["tss"]  # Assign within array.
        res["start"][forward] -= upstream
        res["start"][~forward] -= downstream
        res["start"][res["start"] < 0] = 0
        res["stop"][~forward] += upstream
        res["stop"][forward] += downstream
        res["strand"] = df["strand"]
        res["gene_stable_id"] = df["gene_stable_id"]
        res = pd.DataFrame(res)
        res = res[res["start"] != res["stop"]]  # can happen by the 0 limiting
        return res

    return do_promotorize, [], (upstream, downstream)


def shift(basepairs):
    def do_shift(df):
        res = {
            "chr": df["chr"],
            "start": df["start"] + basepairs,
            "stop": df["stop"] + basepairs,
        }
        return pd.DataFrame(res)

    return do_shift


def summit(summit_annotator):
    def do_summits(df):
        summit_col = summit_annotator.columns[0]
        starts = (df["start"] + df[summit_col]).astype(int)
        res = {"chr": df["chr"], "start": starts, "stop": starts + 1}
        return pd.DataFrame(res)

    return do_summits, [summit_annotator]


def merge_connected():
    """Merge regions that are next to each other.
    100..200, 200..300 becomes 100..300
    """

    def do_merge(df):
        from mbf_nested_intervals import merge_df_intervals

        return merge_df_intervals(df, lambda iv: iv.merge_connected())

    return do_merge


class LiftOver(object):
    def __init__(self):
        import mbf.genomes
        from mbf.externals.kent import LiftOver as LiftOverAlgorithm

        self.data_path = mbf.genomes.data_path / "liftovers"
        self.replacements = {"hg19to38": {"11_gl000202_random": "GL000202.1"}}
        self.algo = LiftOverAlgorithm()

    def do_liftover(self, listOfChromosomeIntervals, chain_file):
        """perform a lift over. Error messages are silently swallowed!"""
        tmp_input = tempfile.NamedTemporaryFile(mode="wb")
        tmp_output = tempfile.NamedTemporaryFile(mode="wb")
        tmp_error = tempfile.NamedTemporaryFile(mode="wb")
        max_len = 0
        listOfChromosomeIntervals = [list(row) for row in listOfChromosomeIntervals]
        for row in listOfChromosomeIntervals:
            tmp_input.write(b" ".join(to_bytes(str(x)) for x in row))
            tmp_input.write(b"\n")
            max_len = max(len(row), max_len)
        tmp_input.write(b"\n")
        tmp_input.flush()  # it's magic ;)
        cmd = [
            self.algo.primary_binary,
            tmp_input.name,
            chain_file,
            tmp_output.name,
            tmp_error.name,
        ]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        dummy_stdout, stderr = p.communicate()
        if p.returncode != 0:  # pragma: no cover
            raise ValueError(
                "do_liftover failed. Returncode: %s, stderr: %s"
                % (p.returncode, stderr)
            )
        tmp_output_in = open(tmp_output.name, "rb")
        res = []
        for row in tmp_output_in:
            row = row.strip().split(b"\t")
            row[0] = to_string(row[0])
            row[1] = int(row[1])
            row[2] = int(row[2])
            res.append(tuple(row))
        tmp_error_in = open(tmp_error.name, "rb")
        tmp_error_in.read()
        tmp_input.close()
        tmp_output.close()
        tmp_error.close()
        return res

    def get_convert_func(self, key, keep_name=False, filter_to_these_chromosomes=None):
        """Note that filter_to_these_chromosomes is after the replacements have kicked in"""
        chain_file = self.data_path / (key + ".over.chain")
        if not chain_file.exists():  # pragma: no cover
            raise ValueError("invalid liftover key, file not found: %s" % chain_file)
        if filter_to_these_chromosomes:
            filter_to_these_chromosomes = set(filter_to_these_chromosomes)

        def do_convert(df):
            if df.index.duplicated().any():  # pragma: no cover
                raise ValueError("liftover only works with unique indices")
            df.index = [str(x) for x in df.index]
            input_tuples = [
                ("chr" + row["chr"], row["start"], row["stop"], idx)
                for idx, row in df.iterrows()
            ]

            output_tuples = self.do_liftover(input_tuples, chain_file)
            output_lists = list(zip(*output_tuples))
            res = pd.DataFrame(
                {
                    "chr": output_lists[0],
                    "start": output_lists[1],
                    "stop": output_lists[2],
                    "parent": [x.decode("utf-8") for x in output_lists[3]],
                }
            ).set_index("parent")
            new_chr = []
            for x in res["chr"]:
                x = x[3:]
                # these are untested as of 2019-03-27
                if x == "m":  # pragma: no cover
                    x = "MT"
                elif (
                    key in self.replacements and x in self.replacements[key]
                ):  # pragma: no cover
                    x = self.replacements[key][x]
                new_chr.append(x)
            res["chr"] = new_chr
            for col in df.columns:
                if col not in res.columns:
                    res = res.assign(**{col: df[col]})
            if filter_to_these_chromosomes:
                res = res[res["chr"].isin(filter_to_these_chromosomes)]
            return res

        if ppg.inside_ppg():
            do_convert.dependencies = [
                ppg.FileTimeInvariant(chain_file),
                ppg.FunctionInvariant(
                    "genomics.regions.convert.LiftOver.do_liftover",
                    LiftOver.do_liftover,
                ),
            ]
        return do_convert


def lift_over(from_to, keep_name=False, filter_to_these_chromosomes=None):
    """Map a genome to another genome.
    from_to looks like hg19ToHg38
    see mbf.genomes/data/liftovers for the list currently supported"""
    return LiftOver().get_convert_func(
        from_to,
        keep_name=keep_name,
        filter_to_these_chromosomes=filter_to_these_chromosomes,
    )


def cookie_cutter(bp):
    """transform all their binding regions to -1/2 * bp ... 1/2 * bp centered
    around the old midpoint... (so pass in the final size of the region)
    inspired by Lupien et al (doi 10.1016/j.cell.2008.01.018")
    """

    def convert(df):
        peak_lengths = df["stop"] - df["start"]
        centers = np.array(df["start"] + peak_lengths // 2, dtype=np.int32)
        new_starts = centers - bp // 2
        new_stops = new_starts + bp
        new_starts[new_starts < 0] = 0
        res = pd.DataFrame({"chr": df["chr"], "start": new_starts, "stop": new_stops})
        if "strand" in df.columns:  # pragma: no branch
            res["strand"] = df["strand"]
        return res

    return convert, [], bp


def cookie_summit(summit_annotator, bp, drop_those_outside_chromosomes=False):
    """transform all their binding regions to -1/2 * bp ... 1/2 * bp centered
    around the summit (so pass in the final size of the region)

    if @drop_those_outside_chromosomes is set, regions < 0 are dropped
    """

    def do_summits(df):
        summit_col = summit_annotator.columns[0]
        res = {
            "chr": df["chr"],
            "start": df["start"] + df[summit_col].astype(int) - bp // 2,
            "stop": df["start"] + df[summit_col].astype(int) + bp // 2,
        }
        res = pd.DataFrame(res)
        if drop_those_outside_chromosomes:
            res = res[res["start"] >= 0]
        else:
            res = res.assign(start=res["start"].clip(lower=0))
        return res

    return do_summits, [summit_annotator], (bp, drop_those_outside_chromosomes)


def windows(window_size, drop_smaller_windows=False):
    """Chuck the region into window_size sized windows.
    if @drop_smaller_windows is True, the right most windows get chopped"""

    def create_windows(df):
        res = {"chr": [], "start": [], "stop": []}
        for dummy_idx, row in df.iterrows():
            for start in range(row["start"], row["stop"], window_size):
                stop = min(start + window_size, row["stop"])
                if drop_smaller_windows and stop - start < window_size:
                    continue
                res["chr"].append(row["chr"])
                res["start"].append(start)
                res["stop"].append(stop)
        return pd.DataFrame(res)

    return create_windows, [], (window_size, drop_smaller_windows)
