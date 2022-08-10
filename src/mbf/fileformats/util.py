from pathlib import Path


def open_file(fileNameOrHandle, mode="rb"):
    """Transparently open compressed or uncompressed files"""
    if hasattr(fileNameOrHandle, "read"):
        return fileNameOrHandle
    elif isinstance(fileNameOrHandle, Path):
        fileNameOrHandle = str(fileNameOrHandle)
    if fileNameOrHandle.endswith(".gz"):
        import gzip

        return gzip.GzipFile(fileNameOrHandle, mode)
    elif fileNameOrHandle.endswith(".bz2"):
        import bz2

        return bz2.BZ2File(fileNameOrHandle, mode)
    else:
        return open(fileNameOrHandle, mode)


def chunkify(handle, separator, block_size=None):
    """take a file handle and split it at separator, reading in efficently in 50 mb blocks or so"""
    if block_size is None:
        block_size = 50 * 1024 * 1024
    chunk = handle.read(block_size)
    chunk = chunk.split(separator)
    while True:
        for k in chunk[:-1]:
            yield k
        next = handle.read(block_size)
        if next:
            chunk = chunk[-1] + next
            chunk = chunk.split(separator)
        else:
            yield chunk[-1]
            break


def pathify(output_filename, default, create_parents=True):
    if output_filename is None:
        res = Path(default)
    else:
        res = Path(output_filename)

    if create_parents:
        res.parent.mkdir(exist_ok=True)
    return res

