"""sequence/chipseq formats"""
from .util import open_file, chunkify
import numpy
import pandas as pd
import os
import tempfile
import subprocess
import re


def normalize_strand(x):
    if isinstance(x, bytes):
        x = x.decode("utf-8")
    if x == 1 or x == 0 or x == -1:
        return x
    if x == "+":
        return 1
    elif x == "-":
        return -1
    else:
        return 0


class BedEntry(object):
    __slots__ = ["refseq", "position", "length", "strand", "score", "name"]

    def __init__(self, chr, chrStart, chrEnd, name=None, strand=None, score=None):
        self.refseq = chr
        pos = int(chrStart)
        self.position = pos
        self.length = int(chrEnd) - pos
        if strand:
            self.strand = normalize_strand(strand)
        else:
            self.strand = normalize_strand(self.length > 0)
        self.score = numpy.nan
        if name is None:
            name = "Noname"
        self.name = name
        if score is not None:
            self.score = score

    def get_read_length(self):
        return self.length

    def __len__(self):
        return self.length

    def __repr__(self):
        return "BedEntry(%s, %i, %i)" % (
            repr(self.refseq),
            self.position,
            self.position + self.length,
        )

    def __str__(self):
        print(
            "Bed Entry chr=%s, start=%i, length=%i, strand=%s, score=%s, name=%s"
            % (
                self.refseq,
                self.position,
                self.length,
                self.strand,
                self.score,
                self.name,
            )
        )


def read_bed(filenameOrFileObject, filter_to_track=None):
    """Read bed file into a list. Either include all tracks (default),
    or filter to specific one with filter_to_track = track_name"""
    res = []
    for e in read_bed_iterator(filenameOrFileObject, filter_to_track):
        res.append(e)
    return res


def read_bed_iterator(filenameOrFileObject, filter_to_track=None):
    fo = open_file(filenameOrFileObject, "rb")
    ever_inside_track = False
    if filter_to_track:
        inside_track = False
    else:
        inside_track = True
        ever_inside_track = True
    for row in chunkify(fo, b"\n"):
        if row.startswith(b"track"):
            if filter_to_track:
                track_name = re.findall(b'name="([^"]+)"', row)[0].decode("utf-8")
                if track_name == filter_to_track:
                    inside_track = True
                    ever_inside_track = True
                else:
                    inside_track = False
            else:
                inside_track = True
                ever_inside_track = True
            # trackInfo = row
            pass
        elif row.startswith(b"#"):  # not really a comment character...
            continue
        else:
            if not inside_track:
                continue
            try:
                if not row:
                    continue
                row = row.split(b"\t")
                e = BedEntry(
                    row[0].decode("utf-8"), row[1], row[2]
                )  # bed does start at 0
                try:
                    e.name = row[3].decode("utf-8")
                    e.score = float(row[4])
                except IndexError:
                    pass
                except ValueError:
                    pass
                try:
                    e.strand = normalize_strand(row[5])
                except IndexError:
                    pass
                yield e
            except Exception as e:
                raise ValueError("Could not parse row: %s" % row, e)
    if not ever_inside_track and filter_to_track:
        raise ValueError("Could not find track %s" % filter_to_track)


def write_bed_header(file_handle, track_name):
    file_handle.write(
        ('track name="%s" description="" useScore=0\n' % track_name).encode("utf-8")
    )


def write_bed_entry_short(file_handle, entry, name_to_chromosome_lookup_or_none=None):
    if name_to_chromosome_lookup_or_none:
        chr = name_to_chromosome_lookup_or_none[entry.refseq]
    else:
        chr = entry.refseq
    out = [
        chr,  # chromosome
        entry.position,  # start
        entry.position + len(entry),  # end
    ]
    file_handle.write("\t".join((str(x) for x in out)) + "\n")


def write_bed_entry_long(file_handle, entry, name_to_chromosome_lookup_or_none=None):
    if name_to_chromosome_lookup_or_none:
        chr = name_to_chromosome_lookup_or_none[entry.refseq]
    else:
        chr = entry.refseq
    out = [
        chr,  # chromosome
        entry.position,  # start
        entry.position + len(entry),  # end
        entry.name,
        "." if numpy.isnan(entry.score) else entry.score,
        "+" if entry.strand == 1 else "-" if entry.strand == -1 else ".",
        entry.position,  # this start
        entry.position + len(entry),
    ]
    file_handle.write(b"\t".join((str(x).encode("utf-8") for x in out)) + b"\n")


def write_bed(
    filenameOrFileObject,
    reads,
    name_to_chromosome_lookup_or_none,
    track_name,
    include_header=True,
    minimal=False,
):
    fo = open_file(filenameOrFileObject, "wb")
    if include_header:
        write_bed_header(fo, track_name)
    if minimal:
        write_method = write_bed_entry_short
    else:
        write_method = write_bed_entry_long
    for read in reads:
        write_method(fo, read, name_to_chromosome_lookup_or_none)
    if fo != filenameOrFileObject:
        fo.close()


def bed_to_bigbed(
    input_bed_filename, output_filename, chromosome_lengths, already_sorted=False
):
    """Convert an existing bed file into bigbed. Chromosome lengths is a dictionary"""
    # todo: actually use the extern bedToBigBed?
    from mbf.externals.kent import BedToBigBed  # noqa: F401

    chrom_sizes_file = tempfile.NamedTemporaryFile(suffix=".sizes")
    for chr, length in sorted(chromosome_lengths.items()):
        chrom_sizes_file.write(("%s\t%i\n" % (chr, length)).encode("utf-8"))
    chrom_sizes_file.flush()

    input = open(input_bed_filename, "rb")
    if not already_sorted:
        missing_chroms = set()
        out = list()
        chrom_lenghts = chromosome_lengths
        for line in input:
            line = line.strip()
            if line and not line.startswith(b"#") and not line.startswith(b"track"):
                line = line.split(b"\t")
                try:
                    chr = line[0].decode("utf-8")
                    start = min(max(0, int(line[1])), chrom_lenghts[chr] - 1)
                    stop = min(int(line[2]), chrom_lenghts[chr] - 1)
                    strand = None
                    name = None
                    try:
                        name = line[3].decode("utf-8")
                        strand = line[5].decode("utf-8")
                    except IndexError:
                        if name is None:
                            name = "."
                        if strand is None:
                            strand = "."
                    out.append(
                        (chr, min(start, stop), max(start, stop), name, 0, strand)
                    )
                except KeyError:
                    missing_chroms.add(chr)
        if missing_chroms:
            print("could not find the following chromosomes %s" % (missing_chroms,))
        out.sort()
        inputtf = tempfile.NamedTemporaryFile()
        for chr, start, stop, name, score, strand in out:
            inputtf.write(
                (
                    "%s\t%s\t%s\t%s\t%s\t%s\n" % (chr, start, stop, name, score, strand)
                ).encode("utf-8")
            )
        inputtf.flush()
        inputtf.seek(0, 0)
        input.close()
    else:
        inputtf = input
    in_name = os.path.abspath(inputtf.name)
    cmd = [
        "bedToBigBed",
        "-tab",
        in_name,
        os.path.abspath(chrom_sizes_file.name),
        os.path.abspath(output_filename),
        "-type=bed6",
    ]
    print(cmd)
    p = subprocess.Popen(cmd)
    p.communicate()
    inputtf.close()
    chrom_sizes_file.close()
    if p.returncode != 0:
        raise ValueError("bedToBigBed returned an error code")


def write_bigbed(
    input_dataframe_or_list_of_bed_entries, output_filename, chromosome_lengths
):
    """Take either a list of BedEntry objects, a dataframe with the slot names of BedEntry objects, or
    a DataFrame with {chr, start, stop, strand, score, name}.
    Either way, strand needs to be one of -1,0,1

    """
    if not isinstance(input_dataframe_or_list_of_bed_entries, pd.DataFrame):
        r = {}
        for column in ["refseq", "position", "length", "strand", "score", "name"]:
            r[column] = []
            for e in input_dataframe_or_list_of_bed_entries:
                r[column].append(getattr(e, column))
        df = pd.DataFrame(r)
    else:
        df = input_dataframe_or_list_of_bed_entries

    if "chr" in df.columns and "start" in df.columns and "stop" in df.columns:
        df = df.sort_values(["chr", "start"], ascending=[True, True])
    elif "refseq" in df.columns and "position" in df.columns and "length" in df.columns:
        df = df.rename(columns={"refseq": "chr", "position": "start"})
        df = df.assign(stop=df["start"] + df["length"])
        df = df.sort_values(["chr", "start"], ascending=[True, True])
    else:
        raise ValueError(
            "This dataframe did not contain the necessary bed columns (either (chr, start, stop) or (refseq, position, length)"
        )
    if "strand" not in df.columns:
        df = df.assign("strand", 0)
    if "name" not in df.columns:
        df = df.assign("name", "Noname")
    if "score" not in df.columns:
        df = df.assign("score", numpy.nan)
    of = tempfile.NamedTemporaryFile(suffix=".bed")
    for dummy_idx, row in df.iterrows():
        output = "%s\t%i\t%i\t%s\t%s\t%s\n" % (
            row["chr"],
            row["start"],
            row["stop"],
            row["name"],
            ("%f" % (row["score"],)) if (row["score"] is not None) else ".",
            ("+" if row["strand"] == 1 else "-" if row["strand"] == -1 else "+"),
        )
        output = output.encode("utf-8")
        of.write(output)
    of.flush()
    bed_to_bigbed(of.name, output_filename, chromosome_lengths)
    of.close()


def read_bigbed(filename, chromosome_lengths, chromosome_mangler=lambda x: x):
    import pyBigWig

    bb = pyBigWig.open(filename)
    chr_lengths = chromosome_lengths
    data = {"chr": [], "start": [], "stop": [], "strand": [], "name": []}
    for chr in chr_lengths:
        cl = chr_lengths[chr]
        chr_to_query = chromosome_mangler(chr)
        if chr_to_query is None:
            continue
        try:
            it = bb.entries(chr_to_query, 0, cl)
        except (
            Exception
        ) as e:  # I have no idea where the RuntimeError exception typecomes from.
            if "Invalid interval bounds" in str(e):
                it = None
            else:
                raise

        if it is None:  # no such chromosome. Tolerable if it's a contig or such.
            # If none of the chromosome names match,
            # we raise later because of an empty big file.
            continue
        for entry in it:
            data["chr"].append(chr)
            data["start"].append(entry[0])
            data["stop"].append(entry[1])
            more = entry[2].split("\t")
            strand = more[2]
            data["strand"].append(1 if strand == "+" else -1 if strand == "-" else 0)
            data["name"].append(more[0])
    bb.close()
    return pd.DataFrame(data)
