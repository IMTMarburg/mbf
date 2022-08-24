# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE=1

# """Collection of fastq-parsers (processors)
# that read fastq (or any other read format)
# filter and massage the reads,
# and write them out again.
# """


import _common
import pypipegraph as ppg
from collections import namedtuple

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans

rev_comp_table = maketrans(
    b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR"
)


def iterate_fastq(fn, reverse_reads):
    op = _common.BlockedFileAdaptor(fn)
    while True:
        try:
            name = op.readline()[1:-1]
            seq = op.readline()[:-1]
            dummy = op.readline()
            qual = op.readline()[:-1]
            if reverse_reads:
                seq = seq[::-1].translate(rev_comp_table)
                qual = qual[::-1]
            yield (seq, qual, name)
        except StopIteration:
            break


def get_iterator(read_creator):
    if read_creator == "fastq":
        return iterate_fastq
    else:
        raise ValueError("Invalid read creator")  # pragma: no cover


class Straight(object):
    """Take a set of fastqs (any old format),
    and prepare an aligner input file.

    This is the base class that neither filters,
    nor truncates reads, it simply creates one big .fastq,
    uncompressed for aligner input.

    See below for more advanced variants"""

    def generate_aligner_input(
        self, output_filename, list_of_fastqs, reverse_reads, read_creator="fastq"
    ):
        our_iter = get_iterator(read_creator)
        with open(output_filename, "wb") as op:
            for fn in list_of_fastqs:
                if reverse_reads or read_creator != "fastq":
                    for seq, qual, name in our_iter(fn, reverse_reads):
                        op.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")
                else:  # no need to reverse the reads, we can just copy/paste
                    for block in _common.read_file_blocked(fn, 20 * 1024 * 1024):
                        op.write(block)

    def get_dependencies(self, output_filename):
        return []


class Filtered(Straight):
    """Filter reads with a callback func that takes seq, 
    qual, name and returns True/False"""

    def __init__(self, filter_func):
        Straight.__init__(self)
        self.filter_func = filter_func

    def generate_aligner_input(
        self, output_filename, list_of_fastqs, reverse_reads, read_creator="fastq"
    ):
        if read_creator == "fastq":
            our_iter = iterate_fastq
        else:
            raise ValueError("Invalid read creator")  # pragma: no cover
        with open(output_filename, "wb") as op:
            for fn in list_of_fastqs:
                for seq, qual, name in our_iter(fn, reverse_reads):
                    if self.filter_func(seq, qual, name):
                        op.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")

    def get_dependencies(self, output_filename):
        return [ppg.FunctionInvariant(output_filename + "_filter", self.filter_func)]


class Paired_Filtered(Straight):
    """Filter reads with a callback func that takes seq1,qual1, name1,
    seq2, qual2, name2 and returns True/False

    All values are bytes!
    """

    def __init__(self, filter_func):
        Straight.__init__(self)
        self.filter_func = filter_func

    def generate_aligner_input_paired(
        self,
        output_filename1,
        output_filename2,
        list_of_fastq_pairs,
        reverse_reads,
        read_creator="fastq",
    ):
        if read_creator == "fastq":
            our_iter = iterate_fastq
        else:
            raise ValueError("Invalid read creator")  # pragma: no cover
        with open(output_filename1, "wb") as op1:
            with open(output_filename2, "wb") as op2:
                for fn1, fn2 in list_of_fastq_pairs:
                    for tup in zip(
                        our_iter(fn1, reverse_reads), our_iter(fn2, reverse_reads)
                    ):
                        seq1, qual1, name1 = tup[0]
                        seq2, qual2, name2 = tup[1]
                        if self.filter_func(seq1, qual1, name1, seq2, qual2, name2):
                            op1.write(
                                b"@" + name1 + b"\n" + seq1 + b"\n+\n" + qual1 + b"\n"
                            )
                            op2.write(
                                b"@" + name2 + b"\n" + seq2 + b"\n+\n" + qual2 + b"\n"
                            )

    def get_dependencies(self, output_filenames):
        return [
            ppg.FunctionInvariant(output_filenames[0] + "_filter", self.filter_func)
        ]


class Paired_Filtered_Trimmed(Straight):
    """Filter reads with a callback func that takes seq1,qual1, name1,
    seq2, qual2, name2 and returns truncated reads/qualities
    """

    def __init__(self, filter_func):
        Straight.__init__(self)
        self.filter_func = filter_func

    def generate_aligner_input_paired(
        self,
        output_filename1,
        output_filename2,
        list_of_fastq_pairs,
        reverse_reads,
        read_creator="fastq",
    ):
        if read_creator == "fastq":
            our_iter = iterate_fastq
        else:
            raise ValueError("Invalid read creator")  # pragma: no cover
        counter = 0
        seen = 0
        with open(output_filename1, "wb") as op1:
            with open(output_filename2, "wb") as op2:
                for fn1, fn2 in list_of_fastq_pairs:
                    for tup in zip(
                        our_iter(fn1, reverse_reads), our_iter(fn2, reverse_reads)
                    ):
                        seq1, qual1, name1 = tup[0]
                        seq2, qual2, name2 = tup[1]
                        seen += 1
                        filtered = self.filter_func(
                            seq1, qual1, name1, seq2, qual2, name2
                        )
                        if filtered is not None:
                            s1, q1, n1, s2, q2, n2 = filtered
                            if s1 is not None and s2 is not None:
                                op1.write(
                                    (b"@" + n1 + b"\n" + s1 + b"\n+\n" + q1 + b"\n")
                                )
                                op2.write(
                                    (b"@" + n2 + b"\n" + s2 + b"\n+\n" + q2 + b"\n")
                                )
                                counter += 1

    def generate_aligner_input(
        self, output_filename, list_of_fastqs, reverse_reads, read_creator="fastq"
    ):
        """This allows to see both mate pairs and select one of them"""

        if read_creator == "fastq":
            our_iter = iterate_fastq
        else:
            raise ValueError("Invalid read creator")  # pragma: no cover
        counter = 0
        seen = 0
        with open(output_filename, "wb") as op:
            for fn1, fn2 in list_of_fastqs:
                for tup in zip(
                    our_iter(fn1, reverse_reads), our_iter(fn2, reverse_reads)
                ):
                    seq1, qual1, name1 = tup[0]
                    seq2, qual2, name2 = tup[1]
                    seen += 1
                    filtered = self.filter_func(
                        seq1, qual1, name1, seq2, qual2, name2
                    )
                    if filtered is not None:
                        s1, q1, n1 = filtered
                        if s1 is not None:
                            op.write(
                                (b"@" + n1 + b"\n" + s1 + b"\n+\n" + q1 + b"\n")
                            )
                            counter += 1

    def get_dependecies(self, output_filenames):
        return [
            ppg.FunctionInvariant(output_filenames[0] + "_filter", self.filter_func)
        ]


class QualityFilter(object):
    """Support for old style quality filters.
    they are functions taking (quality, sequence) (as bytes), and returning
    either true/false (to keep the read as is),
    a positive integer n (keep the first n bases),
    a negative integer n (keep the last |n| bases)
    a tuple (start:stop) - apply python slicing"""

    def __init__(self, quality_filter_func):
        self.filter_func = quality_filter_func

    def generate_aligner_input(
        self, output_filename, list_of_fastqs, reverse_reads, read_creator="fastq"
    ):
        if read_creator == "fastq":
            our_iter = iterate_fastq
        else:
            raise ValueError("Invalid read creator")
        with open(output_filename, "wb") as op:
            for fn in list_of_fastqs:
                for seq, qual, name in our_iter(fn, reverse_reads):
                    r = self.filter_func(qual, seq)
                    if r is False:
                        continue

                    elif isinstance(r, tuple):  # treat as python silce
                        seq = seq[r[0] : r[1]]
                        qual = qual[r[0] : r[1]]
                    elif isinstance(r, int) and r is not True:
                        if r > 0:  # cut right
                            seq = seq[:r]
                            qual = qual[:r]
                        elif r < 0:  # ct left
                            seq = seq[r:]
                            qual = qual[r:]
                        else:
                            raise ValueError(
                                "0 is not a valid qualityfilter return - say True if you want to keep the read as is"
                            )
                    op.write(b"@" + name + b"\n" + seq + b"\n+\n" + qual + b"\n")

    def get_dependencies(self, output_filename):
        return [ppg.FunctionInvariant(output_filename + "_filter", self.filter_func),
                ppg.FunctionInvariant('QualityFilter.generate_aligner_input',QualityFilter.generate_aligner_input)
                ]


def CutAdaptThreePrime(*args, **kwargs):
    raise DeprecationWarning(
        "Use CutAdapt with adapter_sequence_begin=None and adapter_sequence_end=your_seq instead"
    )


def CutAdapt(
    adapter_sequence_begin,
    adapter_sequence_end,
    keep_adapter_less_sequences,
    minimum_remaining_length=10,
    maximal_error_rate=3,
):
    """Search for adapters with cutadapt and trim reads accordingly.
        @adapter_sequence_begin, @adapter_sequence_end may be a sequence, an int (=> python slice,
        positive for begin, negative for end) or None in which case no trimming on that side is performed
        @keep_adapter_less_sequences - if this is False, we throw away reads without adapters
        @minimum_remaining_length - throw away reads shorter than this after trimming
        @maximal_error_rate = allowed mismaches in adapter (integer)
        """
    if isinstance(adapter_sequence_end, int):
        if adapter_sequence_end > 0:
            raise ValueError(
                "adapter_sequence_end needs to be a negative integer, was %s"
                % adapter_sequence_end
            )
    if isinstance(adapter_sequence_begin, int):
        if adapter_sequence_begin > 0:  # pragma: no branch
            raise ValueError(
                "adapter_sequence_begin needs to be a positive integer, was %s"
                % adapter_sequence_begin
            )

    import cutadapt
    import cutadapt.align

    AdapterMatch = namedtuple(
        "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
    )
    where = (
        cutadapt.align.START_WITHIN_SEQ2
        | cutadapt.align.STOP_WITHIN_SEQ2
        | cutadapt.align.STOP_WITHIN_SEQ1
    )
    if isinstance(adapter_sequence_begin, str):
        adapter_begin = cutadapt.align.Aligner(
            adapter_sequence_begin if adapter_sequence_begin else "",
            maximal_error_rate / len(adapter_sequence_begin),
            where,
            wildcard_ref=True,
            wildcard_query=False,
            indel_cost=50000,  # we only want mismatches
        )
    else:
        if adapter_sequence_begin is None:  # pragma: no branch
            adapter_sequence_begin = 0
        adapter_begin = None
    if isinstance(adapter_sequence_end, str):
        adapter_end = cutadapt.align.Aligner(
            adapter_sequence_end if adapter_sequence_end else "",
            maximal_error_rate / len(adapter_sequence_end),
            where,
            wildcard_ref=True,
            wildcard_query=False,
            indel_cost=50000,  # we only want mismatches..
        )
    else:
        adapter_end = None

    def match(adapter, seq):
        alignment = adapter.locate(seq.decode("latin1"))
        if alignment is None:
            return None
        _match = AdapterMatch(*alignment)
        return _match

    def qf(qual, seq):
        if isinstance(adapter_sequence_begin, str):
            match_begin = match(adapter_begin, seq)
            print("begin", match_begin)
            if match_begin is None or match_begin.matches + match_begin.errors < len(
                adapter_sequence_begin
            ):
                if keep_adapter_less_sequences:
                    first_index = 0
                else:
                    return False
            else:
                first_index = match_begin.rstop
        else:
            first_index = adapter_sequence_begin  # slice. None is from start...

        if isinstance(adapter_sequence_end, str):
            match_end = match(adapter_end, seq)
            print("end", match_end)
            if match_end is None or match_end.matches + match_end.errors < len(
                adapter_sequence_end
            ):
                if keep_adapter_less_sequences:
                    second_index = len(seq)
                else:
                    return False
            else:
                second_index = match_end.rstart  # keep till here. positive slice end
        else:
            if adapter_sequence_end is None:
                second_index = len(seq)
            else:
                second_index = adapter_sequence_end
        if second_index < 0:
            second_index = len(seq) + second_index
        print(first_index, second_index)
        if second_index - first_index < minimum_remaining_length:
            return False
        return first_index, second_index

    return QualityFilter(qf)


class UMIExtract(object):
    """Take a set of fastqs 
    and pull out the first N bases as an UMI,
    attach _UMI to the name.
    """

    def __init__(self, umi_length):
        self.umi_length = umi_length

    def generate_aligner_input(
        self, output_filename, list_of_fastqs, reverse_reads, read_creator="fastq"
    ):
        n = self.umi_length
        our_iter = get_iterator(read_creator)
        with open(output_filename, "wb") as op:
            if not reverse_reads:
                for fn in list_of_fastqs:
                    for seq, qual, name in our_iter(fn, False):
                        umi = seq[:n]
                        seq = seq[n:]
                        qual = qual[n:]
                        name = name.split(b" ")
                        name[0] += b"_" + umi
                        op.write(
                            b"@"
                            + b" ".join(name)
                            + b"\n"
                            + seq
                            + b"\n+\n"
                            + qual
                            + b"\n"
                        )
            else:
                raise NotImplementedError("implement for reverse reads")

    def get_dependencies(self, output_filename):
        return [
            ppg.ParameterInvariant(output_filename + "_umiextract", self.umi_length)
        ]


class UMIExtractAndTrim(object):
    """Take a set of fastqs
    and pull out the first N bases as an UMI,
    then trim further @cut_5_prime bases from the front, and
    @cut_3_prime from the end
    attach _UMI to the name.
    """

    def __init__(self, umi_length, cut_5_prime, cut_3_prime=None):
        self.umi_length = umi_length
        self.cut_5_prime = cut_5_prime
        if cut_5_prime < 0:
            raise ValueError("cut_5_prime must be >= 0")
        if cut_3_prime < 0:
            raise ValueError("cut_3_prime must be >= 0")
        if self.cut_5_prime == 0:
            self.cut_5_prime = None
        if cut_3_prime == 0:
            self.cut_3_prime = None
        else:
            self.cut_3_prime = -1 * cut_3_prime

    def generate_aligner_input(
        self, output_filename, list_of_fastqs, reverse_reads, read_creator="fastq"
    ):
        n = self.umi_length
        cut_5 = self.cut_5_prime
        cut_3 = self.cut_3_prime
        our_iter = get_iterator(read_creator)
        with open(output_filename, "wb") as op:
            if not reverse_reads:
                for fn in list_of_fastqs:
                    for seq, qual, name in our_iter(fn, False):
                        umi = seq[:n]
                        seq = seq[n + cut_5 : cut_3]
                        qual = qual[n + cut_5 : cut_3]
                        name = name.split(b" ")
                        name[0] += b"_" + umi
                        op.write(
                            b"@"
                            + b" ".join(name)
                            + b"\n"
                            + seq
                            + b"\n+\n"
                            + qual
                            + b"\n"
                        )
            else:
                raise NotImplementedError("implement for reverse reads")

    def get_dependencies(self, output_filename):
        return [
            ppg.ParameterInvariant(
                output_filename + "_umiextract",
                (self.umi_length, self.cut_5_prime, self.cut_3_prime),
            )
        ]


def QuantSeqFWD():
    """Take a set of Lexogen QuantSeq FWD
    fastqs, extract the UMIs (into the name),
    throw away the next 12b, and the last bp
    """
    return UMIExtractAndTrim(6, 12, 1)


# TODO: KHmer Filter

# TODO: alternative read generating iterators (eland, SRA, fasta, SOAP, bedlane, bednamelane, geraldlane, ELANDLaneIgnoringQualities, Seqlane, MixedFormatLane, fastq igonring qualities  )
