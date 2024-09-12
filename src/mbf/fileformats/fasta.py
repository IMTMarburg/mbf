from .util import chunkify, open_file


def iterate_fasta(filename_or_handle, keyFunc=None):
    """An iterator over a fasta file
    Yields tupples of key, sequence on each iteration
    """
    o = open_file(filename_or_handle)
    try:
        key = ""
        sequence = []
        for chunk in chunkify(o, b"\n>"):
            if chunk:
                key = chunk[: chunk.find(b"\n")].strip()
                if key.startswith(b">"):
                    key = key[1:]
                if keyFunc:
                    key = keyFunc(key)
                if chunk.find(b"\n") != -1:
                    seq = (
                        chunk[chunk.find(b"\n") + 1 :]
                        .replace(b"\r", b"")
                        .replace(b"\n", b"")
                    )
                else:
                    seq = ""
                yield (key.decode("utf-8"), seq.decode("utf-8"))
        return

        for line in o:
            if line.startswith(b">"):
                if key != "" and len(sequence) > 0:
                    yield (
                        key,
                        b"".join(sequence).replace(b"\n", b"").replace(b"\r", b""),
                    )
                key = line[1:].strip()
                sequence = []
            else:
                sequence.append(line)
        if key != "" and len(sequence) > 0:
            yield (
                key.decode("utf-8"),
                (b"".join(sequence).replace(b"\n", b"").replace(b"\r", b"")).decode(
                    "utf-8"
                ),
            )
    finally:
        o.close()




def wrappedIterator(width):
    def inner(text):
        i = 0
        length = len(text)
        while i < length:
            yield text[i: i + width]
            i += width
    return inner


def dict_to_fasta(fastaDict, filename, doWrap=80, doUpper=True):
    """Writes a Fasta file from a dictionary of
    keys: sequences. If the filename ends with '.gz',
    the output is gzipped.
    Wraps if doWrap is set (at position 80)"""
    genToFasta(fastaDict.items(), filename, doWrap, doUpper)


def gen_to_fasta(gen, filename, doWrap=80, doUpper=True):
    """Take a generator creating (key, sequence) tuples,
    write it to a file. Wraps if doWrap is set.
    @doUpper may be True, then call .upper, it may be False
    then call .lower, or it may be anything else - then keep them as they are.
    """
    o = open_file(filename, 'wb')
    for key, value in gen:
        if hasattr(key, 'encode'):
            key = key.encode('utf-8')
        if hasattr(value, 'encode'):
            value = value.encode('utf-8')
        o.write(b'>' + key + b"\n")
        if doWrap:
            it = wrappedIterator(doWrap)
            for line in it(value):
                if doUpper:
                    o.write(line.upper() + b"\n")
                elif doUpper is False:
                    o.write(line.lower() + b"\n")
                else:
                    o.write(line + b"\n")
        else:
            o.write(value + b"\n")
    o.close()


