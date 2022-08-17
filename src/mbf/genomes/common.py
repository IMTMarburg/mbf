from mbf.fileformats.util import open_file, chunkify


def iter_fasta(filenameOrFileLikeObject, keyFunc=None, block_size=None):
    """An iterator over a fasta file (raw or gzipped).
    Yields tupples of key, sequence  (bytes!) on each iteration
    """
    o = open_file(filenameOrFileLikeObject)
    key = ""
    for chunk in chunkify(o, b"\n>", block_size=block_size):
        key = chunk[: chunk.find(b"\n")].strip()
        if key.startswith(b">"):
            key = key[1:]
        if keyFunc:
            key = keyFunc(key)
        if chunk.find(b"\n") != -1:
            seq = chunk[chunk.find(b"\n") + 1 :].replace(b"\r", b"").replace(b"\n", b"")
        else:
            raise ValueError("Should not be reached")  # pragma: no cover
            # seq = b""
        yield (key, seq)
    return


def wrappedIterator(width):
    def inner(text):
        i = 0
        length = len(text)
        while i < length:
            yield text[i : i + width]
            i += width

    return inner


rc_table = str.maketrans("agctAGCT", "tcgaTCGA")
iupac_forward = "ACGTRYMKSWBDHVN"
iupac_reverse = "TGCAYRKMSWVHDBN"
iupac_forward += iupac_forward.lower()
iupac_reverse += iupac_reverse.lower()
iupac_rc_table = str.maketrans(
    iupac_forward + iupac_forward.upper(), iupac_reverse + iupac_reverse.upper()
)


def reverse_complement(s):
    """return complementary, reversed sequence to x (keeping case)"""
    return s.translate(rc_table)[::-1]


def reverse_complement_iupac(s):
    """return complementary, reversed sequence to x (keeping case)"""
    return s.translate(iupac_rc_table)[::-1]


universal_genenetic_code = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}


class GeneticCode:
    @classmethod
    def translate_dna(cls, sequence, raise_on_non_multiple_of_three=True):
        if raise_on_non_multiple_of_three and len(sequence) % 3 != 0:
            raise ValueError("len(sequence) was not a multiple of 3")
        genetic_code = cls.genetic_code
        proteinseq = ""
        sequence = sequence.upper()
        if sequence[:3] in cls.start_codons:
            proteinseq += "M"
        else:
            proteinseq += genetic_code[sequence[:3]]
        for n in range(3, len(sequence), 3):
            proteinseq += genetic_code[sequence[n : n + 3]]
        return proteinseq

    @classmethod
    def translate_dna_till_stop(cls, sequence, genetic_code=None):
        genetic_code = cls.genetic_code
        proteinseq = ""
        sequence = sequence.upper()
        sequence = sequence.upper()
        if sequence[:3] in cls.start_codons:
            proteinseq += "M"
        else:
            proteinseq += genetic_code[sequence[:3]]
        for n in range(3, len(sequence), 3):  # pragma: no branch
            try:
                codon = sequence[n : n + 3]
                x = genetic_code[codon]
                proteinseq += x
                if x == "*":
                    break
            except KeyError:
                if len(codon) < 3:
                    raise ValueError("No stop codon found")
                else:
                    raise NotImplementedError(  # pragma: no cover
                        "Incomplete genetic code?, codon %s not found" % codon
                    )
        return proteinseq


class EukaryoticCode(GeneticCode):
    """Genetic code for eukaryotes"""

    genetic_code = universal_genenetic_code
    start_codons = ["ATG"]


class ProkaryoticCode(GeneticCode):
    """Genetic code for prokaryotes - e.g. E. coli"""

    genetic_code = universal_genenetic_code
    # for e coli, from wikipedia) - 83% AUG (3542/4284), 14% (612) GUG, 3% (103) UUG[7] and one or two others (e.g., an AUU and possibly a CUG
    start_codons = ["ATG", "GTG", "TTG", " ATT", "CTG"]


def df_to_rows(df, columns_to_choose=None):
    """Turn a DataFrame into named tuples
    index -> {columnA: X, columnB: Y}
    You can then use that much faster than
    accessing df.loc[] repeatedly
    right now (2019-03-25) for whatever reason
    """
    if columns_to_choose is None:  # pragma: no branch
        pass  # pragma: no cover
    else:
        df = df[columns_to_choose]
    if df.index.duplicated().any():  # pragma: no cover
        raise ValueError("df_to_rows needs a unique index")
    result = {}
    for row in df.itertuples():
        result[row[0]] = row
    return result
