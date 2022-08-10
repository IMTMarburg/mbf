"""Tests for the various combinations you can evoke
our cython based fastq reader with."""
import os
import pytest
import pypipegraph as ppg
from pathlib import Path

from mbf_align import fastq2
from mbf_sampledata import get_sample_data


@pytest.mark.usefixtures("new_pipegraph")
class TestFastqS:
    def test_straight_copy(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_2
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"""
        try:
            with open("test_straight_copy.fastq", "wb") as op:
                op.write(test)
            x = fastq2.Straight()
            x.generate_aligner_input(
                "test_straight_copy.out.fastq", ["test_straight_copy.fastq"], False
            )
            with open("test_straight_copy.out.fastq", "rb") as op:
                was = op.read()
                assert test == was

        finally:
            if os.path.exists("test_straight_copy.out.fastq"):
                os.unlink("test_straight_copy.out.fastq")

    def test_straight_reverse(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_2
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCA
"""
        should = b"""@SEQ_ID_1
AAACTGTGAGTTGAACAAATGGATTTACTATTTGATCGATACTGCTTTGAACCCCAAATT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAcA
@SEQ_ID_2
AAACTGTGAGTTGAACAAATGGATTTACTATTTGATCGATACTGCTTTGAACCCCAAATA
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
AAACTGTGAGTTGAACAAATGGATTTACTATTTGATCGATACTGCTTTGAACCCCATTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_2
AAACTGTGAGTTGAACAAATGGATTTACTATTTGATCGATACTGCTTTGAACCCCAAATG
+
ACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"""
        try:
            with open("test_straight_reverse.fastq", "wb") as op:
                op.write(test)
            x = fastq2.Straight()
            x.generate_aligner_input(
                "test_straight_reverse.out.fastq", ["test_straight_reverse.fastq"], True
            )
            with open("test_straight_reverse.out.fastq", "rb") as op:
                was = op.read()
                assert should == was

        finally:
            if os.path.exists("test_straight_reverse.out.fastq"):
                os.unlink("test_straight_reverse.out.fastq")

    def test_filter(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_4
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCA
"""
        should = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""

        try:
            with open("test_filter.fastq", "wb") as op:
                op.write(test)

            def f(seq, qual, name):
                print(seq, qual, name)
                result = seq.startswith(b"A")
                return result

            x = fastq2.Filtered(f)
            x.generate_aligner_input(
                "test_filter.out.fastq", ["test_filter.fastq"], False
            )
            with open("test_filter.out.fastq", "rb") as op:
                was = op.read()
                assert should == was

        finally:
            if os.path.exists("test_filter.out.fastq"):
                os.unlink("test_filter.out.fastq")

    def test_quality_filter_bool(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_4
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCA
"""
        should = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""

        try:
            with open("test_quality_filter_bool.fastq", "wb") as op:
                op.write(test)

            def f(qual, seq):
                return seq.startswith(b"A")

            x = fastq2.QualityFilter(f)
            x.generate_aligner_input(
                "test_quality_filter_bool.out.fastq",
                ["test_quality_filter_bool.fastq"],
                False,
            )
            with open("test_quality_filter_bool.out.fastq", "rb") as op:
                was = op.read()
                assert should == was

        finally:
            if os.path.exists("test_quality_filter_bool.out.fastq"):
                os.unlink("test_quality_filter_bool.out.fastq")

    def test_quality_filter_positive_integer(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_4
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCA
"""
        should = b"""@SEQ_ID_1
AATTT
+
AcAAA
@SEQ_ID_3
AAAAT
+
BBBBB
"""

        try:
            with open("test_quality_filter_positive_integer.fastq", "wb") as op:
                op.write(test)

            def f(qual, seq):
                if seq.startswith(b"A"):  # keep the first 5 bases
                    return 5
                else:
                    return False

            x = fastq2.QualityFilter(f)
            x.generate_aligner_input(
                "test_quality_filter_positive_integer.out.fastq",
                ["test_quality_filter_positive_integer.fastq"],
                False,
            )
            with open("test_quality_filter_positive_integer.out.fastq", "rb") as op:
                was = op.read()
                assert should == was

        finally:
            if os.path.exists("test_quality_filter_positive_integer.out.fastq"):
                os.unlink("test_quality_filter_positive_integer.out.fastq")

    def test_quality_filter_bool_cut_3(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_4
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCA
"""
        should = b"""@SEQ_ID_1
AGTTT
+
AAAAA
@SEQ_ID_3
AGTTT
+
BBBBB
"""

        try:
            with open("test_quality_filter_bool_cut_3.fastq", "wb") as op:
                op.write(test)

            def f(qual, seq):
                if seq.startswith(b"A"):
                    return -5  # keep the last |-5| bases
                else:
                    return False

            x = fastq2.QualityFilter(f)
            x.generate_aligner_input(
                "test_quality_filter_bool_cut_3.out.fastq",
                ["test_quality_filter_bool_cut_3.fastq"],
                False,
            )
            with open("test_quality_filter_bool_cut_3.out.fastq", "rb") as op:
                was = op.read()
                assert should == was

        finally:
            if os.path.exists("test_quality_filter_bool_cut_3.out.fastq"):
                os.unlink("test_quality_filter_bool_cut_3.out.fastq")

    def test_quality_filter_tuple(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBcBBBBB
@SEQ_ID_4
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCA
"""
        should = b"""@SEQ_ID_1
GGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCAC
+
ghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
GGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCAC
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBc
"""

        try:
            with open("test_quality_filter_tuple.fastq", "wb") as op:
                op.write(test)

            def f(qual, seq):
                if seq.startswith(b"A"):
                    return (6, -5)
                else:
                    return False

            x = fastq2.QualityFilter(f)
            x.generate_aligner_input(
                "test_quality_filter_tuple.out.fastq",
                ["test_quality_filter_tuple.fastq"],
                False,
            )
            with open("test_quality_filter_tuple.out.fastq", "rb") as op:
                was = op.read()
                assert should == was

        finally:
            if os.path.exists("test_quality_filter_tuple.out.fastq"):
                os.unlink("test_quality_filter_tuple.out.fastq")

    def test_cutadapt(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TACCCGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""
        should = b"""@SEQ_ID_1
TTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

        try:
            with open("test_cutadapt.fastq", "wb") as op:
                op.write(test)

            x = fastq2.CutAdapt("AATTTGGGG", None, False, maximal_error_rate=1)
            x.generate_aligner_input(
                "test_cutadapt.out.fastq", ["test_cutadapt.fastq"], False
            )
            with open("test_cutadapt.out.fastq", "rb") as op:
                was = op.read()
                print("should")
                print(should)
                print("was")
                print(was)
                assert should == was

        finally:
            if os.path.exists("test_cutadapt.out.fastq"):
                os.unlink("test_cutadapt.out.fastq")

    def test_cutadapt_both_ends(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TACCCGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""
        should = b"""@SEQ_ID_1
TTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAAC
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

        try:
            with open("test_cutadapt.fastq", "wb") as op:
                op.write(test)

            x = fastq2.CutAdapt("AATTTGGGG", "TCAC", False, maximal_error_rate=1)
            x.generate_aligner_input(
                "test_cutadapt.out.fastq", ["test_cutadapt.fastq"], False
            )
            with open("test_cutadapt.out.fastq", "rb") as op:
                was = op.read()
                print("should")
                print(should)
                print("was")
                print(was)
                assert should == was

        finally:
            if os.path.exists("test_cutadapt.out.fastq"):
                os.unlink("test_cutadapt.out.fastq")

    def test_cutadapt_both_ends_keep(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACAAAAAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""
        should = b"""@SEQ_ID_1
TTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAAC
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACAAAAAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

        try:
            with open("test_cutadapt.fastq", "wb") as op:
                op.write(test)

            x = fastq2.CutAdapt("AATTTGGGG", "TCACAGTTT", True, maximal_error_rate=1)
            x.generate_aligner_input(
                "test_cutadapt.out.fastq", ["test_cutadapt.fastq"], False
            )
            with open("test_cutadapt.out.fastq", "rb") as op:
                was = op.read()
                print("should")
                print(should)
                print("was")
                print(was)
                assert should == was

        finally:
            if os.path.exists("test_cutadapt.out.fastq"):
                os.unlink("test_cutadapt.out.fastq")

    def test_cutadapt_both_ends_no_keep(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_3
AATTTGGGGAAAATCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACAAAAAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""
        should = b"""@SEQ_ID_1
TTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAAC
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

        try:
            with open("test_cutadapt.fastq", "wb") as op:
                op.write(test)

            x = fastq2.CutAdapt("AATTTGGGG", "TCACAGTTT", False, maximal_error_rate=1)
            x.generate_aligner_input(
                "test_cutadapt.out.fastq", ["test_cutadapt.fastq"], False
            )
            with open("test_cutadapt.out.fastq", "rb") as op:
                was = op.read()
                print("should")
                print(should)
                print("was")
                print(was)
                assert should == was

        finally:
            if os.path.exists("test_cutadapt.out.fastq"):
                os.unlink("test_cutadapt.out.fastq")

    def test_cutadapt_only_end(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TACCCGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""
        should = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAAC
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TACCCGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAAC
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

        try:
            with open("test_cutadapt.fastq", "wb") as op:
                op.write(test)

            x = fastq2.CutAdapt(None, "TCAC", False, maximal_error_rate=1)
            x.generate_aligner_input(
                "test_cutadapt.out.fastq", ["test_cutadapt.fastq"], False
            )
            with open("test_cutadapt.out.fastq", "rb") as op:
                was = op.read()
                print("should")
                print(should)
                print("was")
                print(was)
                assert should == was

        finally:
            if os.path.exists("test_cutadapt.out.fastq"):
                os.unlink("test_cutadapt.out.fastq")

    def test_cutadapt_mismatch_straight_keep(self):
        test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AcAAAAghAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TACCCGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""
        should = b"""@SEQ_ID_1
TTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAG
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TACCCGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAG
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

        try:
            with open("test_cutadapt.fastq", "wb") as op:
                op.write(test)

            x = fastq2.CutAdapt("AATTTGGGT", -3, True, maximal_error_rate=1)
            x.generate_aligner_input(
                "test_cutadapt.out.fastq", ["test_cutadapt.fastq"], False
            )
            with open("test_cutadapt.out.fastq", "rb") as op:
                was = op.read()
                print("should")
                print(should)
                print("was")
                print(was)
                assert should == was

        finally:
            if os.path.exists("test_cutadapt.out.fastq"):
                os.unlink("test_cutadapt.out.fastq")


def test_filtered_paired(new_pipegraph):
    import gzip

    r1_name_found = [False]
    r1_qual_found = [False]
    r2_name_found = [False]
    r2_qual_found = [False]

    def f(seq1, qual1, name1, seq2, qual2, name2):
        if name1 == b"HWI-C00113:209:HJCNTBCX2:2:2206:9418:13942 1:N:0:ATCACG":
            r1_name_found[0] = True
        if qual1 == b"DDDDDIIIIIIIIIIIIIIIIIIIIII/<FHHIII<<CHHGHHIHHIIIIH":
            r1_qual_found[0] = True
        if name2 == b"HWI-C00113:209:HJCNTBCX2:2:2206:10802:17968 1:N:0:ATCACG":
            r2_name_found[0] = True
        if qual2 == b"DDDDDIIIIIIIIIIIIIIII1<FHHI/<GHIIII/<EHIIHII/<DHIID":
            r2_qual_found[0] = True
        return seq1.startswith(b"G") and seq2.startswith(b"G")

    x = fastq2.Paired_Filtered(f)
    of1 = "output_R1.fastq"
    of2 = "output_R2.fastq"
    tf1 = open("input_R1_.fastq", "wb")
    tf2 = open("input_R2_.fastq", "wb")
    with gzip.GzipFile(
        get_sample_data(Path("mbf_align/sample_b") / "a_R1_.fastq.gz"), "rb"
    ) as op_in:
        tf1.write(op_in.read())
        tf1.flush()
    with gzip.GzipFile(
        get_sample_data(Path("mbf_align/sample_b") / "a_R2_.fastq.gz"), "rb"
    ) as op_in:
        tf2.write(op_in.read())
        tf2.flush()

    x.generate_aligner_input_paired(of1, of2, [(tf1.name, tf2.name)], False)
    assert r1_name_found[0]
    assert r2_name_found[0]
    assert r1_qual_found[0]
    assert r2_qual_found[0]
    actual1 = Path(of1).read_text()
    actual2 = Path(of2).read_text()
    assert actual1.count("\n") == 4
    assert "@HWI-C00113:209:HJCNTBCX2:2:2206:9559:13855 1:N:0:ATCACG" in actual1
    assert "GCCCAATGTTCGAAATTGCTATTCTACGACAAGGTGCCAGATCTCATCTGA" in actual1
    assert actual2.count("\n") == 4
    assert "@HWI-C00113:209:HJCNTBCX2:2:2206:11052:17798 1:N:0:ATCACG" in actual2
    assert "GTCGGTCCTGAGAGATGGGCGGGCGCCGTTCCGAAAGTACGGGCGATGGCC" in actual2


def test_filtered_depends_on_function_invariant(new_pipegraph):
    def f(seq1, qual1, name1, seq2, qual2, name2):
        return True

    x = fastq2.Paired_Filtered(f)
    deps = x.get_dependencies(["test_R1_.fastq", "test_R2_.fastq"])
    assert isinstance(deps[0], ppg.FunctionInvariant)
    assert "test_R1_.fastq" in deps[0].job_id
    assert deps[0].function is f


def test_filtered_paired_depends_on_function_invariant(new_pipegraph):
    def f(name, seq, qual):
        return True

    x = fastq2.Filtered(f)
    deps = x.get_dependencies("test.fastq")
    assert isinstance(deps[0], ppg.FunctionInvariant)
    assert "test.fastq" in deps[0].job_id
    assert deps[0].function is f


def test_quality_filter_depends_on_function_invariant(new_pipegraph):
    def f(qual, seq):
        return True

    x = fastq2.QualityFilter(f)
    deps = x.get_dependencies("test.fastq")
    assert isinstance(deps[0], ppg.FunctionInvariant)
    assert deps[0].function is f


def test_read_creator_must_be_fastq_right_now(new_pipegraph):
    with pytest.raises(ValueError):
        fastq2.Straight().generate_aligner_input(
            "test.fastq",
            [str(get_sample_data(Path("mbf_align/sample_a)/a.fastq")))],
            False,
            "fail",
        )
    with pytest.raises(ValueError):
        fastq2.Filtered(lambda seq, qual, name: True).generate_aligner_input(
            "test.fastq",
            [str(get_sample_data(Path("mbf_align/sample_a)/a.fastq")))],
            False,
            "fail",
        )
    with pytest.raises(ValueError):
        fastq2.QualityFilter(lambda qual, seq: True).generate_aligner_input(
            "test.fastq",
            [str(get_sample_data(Path("mbf_align/sample_a)/a.fastq")))],
            False,
            "fail",
        )


def test_quality_raises_on_0_return(new_pipegraph):
    with pytest.raises(ValueError):
        fastq2.QualityFilter(lambda qual, seq: 0).generate_aligner_input(
            "test.fastq",
            [str(get_sample_data(Path("mbf_align/sample_a/a.fastq")))],
            False,
        )


def test_cutadapt_three_prime_deprecated():
    with pytest.raises(DeprecationWarning):
        fastq2.CutAdaptThreePrime()


def test_cutadapt_raises_on_negative_adapter_sequences():
    with pytest.raises(ValueError):
        fastq2.CutAdapt(-1, 5, True)
    with pytest.raises(ValueError):
        fastq2.CutAdapt(1, -5, True)


def test_umi_extract(new_pipegraph):
    test = b"""@SEQ_ID_1 123
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_2
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"""
    should = b"""@SEQ_ID_1_AAT 123
TTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2_TAT
TTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2_AAA
ATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_2_CAT
TTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"""

    try:
        with open("test_straight_copy.fastq", "wb") as op:
            op.write(test)
        x = fastq2.UMIExtract(3)
        assert len(x.get_dependencies("test_straight_copy.out.fastq")) == 1
        x.generate_aligner_input(
            "test_straight_copy.out.fastq", ["test_straight_copy.fastq"], False
        )
        with open("test_straight_copy.out.fastq", "rb") as op:
            was = op.read()
            assert should == was

    finally:
        if os.path.exists("test_straight_copy.out.fastq"):
            os.unlink("test_straight_copy.out.fastq")


def test_quantseq_fwd(new_pipegraph):
    test = b"""@SEQ_ID_1
AATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
@SEQ_ID_2 1:N:0:7
TATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2
AAAATGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_2
CATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"""
    should = b"""@SEQ_ID_1_AATTTG
GTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTT
+
GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
@SEQ_ID_2_TATTTG 1:N:0:7
GTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTT
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@SEQ_ID_2_AAAATG
GTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTT
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@SEQ_ID_2_CATTTG
GTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
"""

    try:
        with open("test_straight_copy.fastq", "wb") as op:
            op.write(test)
        x = fastq2.QuantSeqFWD()
        assert len(x.get_dependencies("test_straight_copy.out.fastq")) == 1
        x.generate_aligner_input(
            "test_straight_copy.out.fastq", ["test_straight_copy.fastq"], False
        )
        with open("test_straight_copy.out.fastq", "rb") as op:
            was = op.read()
            assert should == was

    finally:
        if os.path.exists("test_straight_copy.out.fastq"):
            os.unlink("test_straight_copy.out.fastq")
