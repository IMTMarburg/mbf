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
