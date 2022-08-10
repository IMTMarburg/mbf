import numpy as np
import os
import pandas as pd
import io

from .util import open_file


def wiggle_to_intervals(
    filenameOrFileObject, comment_char=None, chromosome_mangler=None
):
    """read a (variable/fixed step) wiggle file and turn it into a df of {chr, start,stop, score}.
    """
    fo = open_file(filenameOrFileObject)
    d = fo.read()
    modes = set()
    if b"variableStep" in d:
        modes.add("variableStep")
    if b"fixedStep" in d:
        modes.add("fixedStep")
    if b"Bed mode" in d:
        modes.add("bed")
    if len(modes) > 1:
        fo.seek(0, os.SEEK_SET)
        return wiggle_to_intervals_slow(
            filenameOrFileObject, comment_char, chromosome_mangler
        )
    lines = np.array(d.strip().split(b"\n"))
    del d
    if comment_char:
        lines = lines[
            ~np.char.startswith(lines, "track")
            & ~np.char.startswith(lines, comment_char)
        ]
    else:
        lines = lines[~np.char.startswith(lines, b"track")]
    if len(modes) == 0:
        if len(lines[0].split()) == 4:
            modes.add("bed")
    if len(modes) == 0:
        raise ValueError(
            "Did not know how to handle %s, no mode found" % filenameOrFileObject
        )
    mode = list(modes)[0]
    if chromosome_mangler is None:
        chromosome_mangler = lambda x: x  # NOQA
    if mode == "bed":
        s = io.BytesIO(b"\n".join(lines))
        recarray = np.loadtxt(
            s,
            dtype=[
                ("chr", "|S50"),
                ("start", np.int32),
                ("stop", np.uint32),
                ("score", np.float),
            ],
        )
        res = pd.DataFrame(
            {
                "chr": [chromosome_mangler(x) for x in recarray["chr"]],
                "start": recarray["start"],
                "stop": recarray["stop"],
                "score": recarray["score"],
            }
        )
        return res
    elif mode == "variableStep":
        borders = np.char.startswith(lines, "variableStep")
        # n umpy.array([x.startswith('variableStep') for x in lines],
        # dtype=np.bool)
        border_offsets = np.where(borders)[0]
        dfs_by_chr = []
        for ii in range(0, len(border_offsets) - 1):
            start = border_offsets[ii]
            stop = border_offsets[ii + 1]
            data_lines = lines[start + 1 : stop]
            row = lines[start]
            chr = row[row.find(b"chrom=") + len(b"chrom=") :]
            chr = chr[: chr.find(" ")]
            chr = chromosome_mangler(chr)
            if row.find(b"span") != -1:
                span = int(row[row.find(b"span=") + len(b"span=") :])
            else:
                span = 1
            just_one_row = len(data_lines) == 1
            data_stream = io.StringIO("\n".join(data_lines) + "\n")
            data_lines = np.loadtxt(data_stream)
            if just_one_row:
                data_lines = np.array([data_lines])
            try:
                starts = np.array(data_lines[:, 0], dtype=np.int)
            except IndexError:
                print(row)
                print(data_lines)
                print(data_lines.value)
                raise
            stops = starts + span
            scores = data_lines[:, 1]
            chrs = np.array([chr] * len(scores), dtype=np.object)
            dfs_by_chr.append(
                pd.DataFrame(
                    {"chr": chrs, "start": starts, "stop": stops, "score": scores}
                )
            )
        # .sort_by(('chr','start'), [True, True]) # we rely on the wiggle being sorte.d..
        return pd.concat(dfs_by_chr)
    elif mode == "fixedStep":
        borders = np.char.startswith(lines, "fixedStep")
        # borders = np.array([x.startswith('fixedStep') for x in lines], dtype=np.bool)
        border_offsets = np.where(borders)[0]
        dfs_by_chr = []
        for ii in range(0, len(border_offsets) - 1):
            data_lines = lines[start + 1 : stop]
            # parse header line...
            row = lines[start]
            chr = row[row.find(b"chrom=") + len(b"chrom=") :]
            chr = chr[: chr.find(" ")]
            chr = chromosome_mangler(chr)
            if row.find(b"span") != -1:
                span = int(row[row.find(b"span=") + len(b"span=") :])
            else:
                raise ValueError("No span in fixed step?!")
            if row.find(b"start") != -1:
                start = row[row.find(b"start=") + len(b"start=") :]
                start = int(start[: start.find(" ")])
            else:
                raise ValueError("No start in fixed step?!")
            if row.find(b"step") != -1:
                step = row[row.find(b"step=") + len(b"step=") :]
                step = int(step[: step.find(" ")])
            else:
                raise ValueError("No step in fixed step?!")
            if step != span:
                raise ValueError(
                    "Parser currently only supports step == span, and they were not equal"
                )
            data_lines = io.StringIO("\n".join(data_lines))
            data_lines = np.loadtxt(data_lines)
            scores = data_lines
            starts = np.array([start + step * ii for ii in range(0, len(scores))])
            stops = stops + span
            chrs = [chr] * len(scores)
            dfs_by_chr.append(
                pd.DataFrame(
                    {"chr": chrs, "start": starts, "stop": stops, "score": scores}
                )
            )
        return pd.concat(dfs_by_chr)
    else:
        raise NotImplementedError("Unknown Wiggle mode")


def wiggle_to_intervals_slow(
    filenameOrFileObject, comment_char=None, chromosome_mangler=None
):
    fo = open_file(filenameOrFileObject)
    result = {}
    mode = None
    chr = None
    span = None
    start = None
    stop = None
    lastScore = None
    data = {"start": [], "stop": [], "score": [], "chr": []}
    if chromosome_mangler is None:
        chromosome_mangler = lambda x: x  # NOQA
    bytes_so_far = 0
    for row in fo:
        bytes_so_far += len(row)
        if (comment_char and row.startswith(comment_char)) or row.startswith("track"):
            continue
        elif row.startswith("variableStep"):
            mode = "variable"
            chr = row[row.find("chrom=") + len("chrom=") :]
            chr = chr[: chr.find(" ")]
            if chr not in result:
                result[chr] = []
            if row.find("span") != -1:
                span = row[row.find("span=") + len("span=") :]
                span = int(span)
            else:
                span = 1
            start = None
            stop = None
            lastScore = None
            continue
        elif row.startswith("fixedStep"):
            mode = "fixed"
            chr = row[row.find("chrom=") + len("chrom=") :]
            chr = chr[: chr.find(" ")]
            if chr not in result:
                result[chr] = []
            if row.find("span") != -1:
                span = row[row.find("span=") + len("span=") :]
                span = int(span)
            else:
                raise ValueError("No span in fixed step?!")
            if row.find("start") != -1:
                start = row[row.find("start=") + len("start=") :]
                start = start[: start.find(" ")]
                start = int(start)
            else:
                raise ValueError("No start in fixed step?!")
            if row.find("step") != -1:
                step = row[row.find("step=") + len("step=") :]
                step = step[: step.find(" ")]
                step = int(step)
            else:
                raise ValueError("No step in fixed step?!")
            if step != span:
                raise ValueError(
                    "Parser currently only supports step == span, and they were not equal"
                )
            lastScore = None
            stop = start
            continue

        elif "Bed format" in row:
            mode = "bed"
            return _wiggle_with_bed_to_intervals(fo, "\t")
        if (
            not mode
        ):  # no mode set so far, perhaps set one now. Not 'continue', we need this row!
            if len(row.split("\t")) == 4:
                mode = "defined"
                # TODO: replace with dataframe csv reading in this case, should
                # be a bit faster
                # because the iterator reads in blocks...
                fo.seek(bytes_so_far - len(row), os.SEEK_SET)
                df = _wiggle_with_bed_to_intervals(fo, "\t")
                df.convert_type("chr", chromosome_mangler)
                return df
            elif len(row.split()) == 4:
                fo.seek(bytes_so_far - len(row), os.SEEK_SET)
                df = _wiggle_with_bed_to_intervals(fo, None)
                df.convert_type("chr", chromosome_mangler)
                return df
            else:
                raise ValueError(
                    "no step mode definied in file. add variablestep/fixedstep to %s"
                    % fo.name
                )
        if mode == "variable":
            row = row.split("\t")
            pos = int(row[0])
            score = float(row[1])
            if score != lastScore:
                if start is not None:
                    data["chr"].append(chromosome_mangler(chr))
                    data["start"].append(start)
                    data["stop"].append(stop)
                    data["score"].append(lastScore)
                lastScore = score
                start = pos
                stop = pos + span
            else:
                stop = pos + span
        elif mode == "fixed":
            score = float(row.strip())
            if score != lastScore:
                if lastScore:
                    data["chr"].append(chromosome_mangler(chr))
                    data["start"].append(start)
                    data["stop"].append(stop)
                    data["score"].append(lastScore)
                start = stop + 1
                lastScore = score
            stop += step
        elif mode == "defined":
            row = row.split("\t")
            chr = row[0]
            start = int(row[1])
            stop = int(row[2])
            score = float(row[3])
            data["chr"].append(chr)
            data["start"].append(start)
            data["stop"].append(stop)
            data["score"].append(score)
    if mode == "variable" or mode == "fixed":
        data["chr"].append(chromosome_mangler(chr))
        data["start"].append(start)
        data["stop"].append(stop)
        data["score"].append(lastScore)
    return pd.DataFrame(data)


def _wiggle_with_bed_to_intervals(file_handle, separator):
    print("calling _wiggle_with_bed_to_intervals")
    res = pd.read_csv.read(file_handle, header=None, sep=separator)
    res.rename_column("column0", "chr")
    res.rename_column("column1", "start")
    res.rename_column("column2", "stop")
    res.rename_column("column3", "score")
    return res

    pass
