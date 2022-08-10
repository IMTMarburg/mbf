import pandas as pd


def read_pandas(filename):
    filename = str(filename)
    if filename.endswith(".xls") or filename.endswith(".xlsx"):
        from xlrd import XLRDError

        try:
            filein = pd.read_excel(filename)
        except XLRDError:
            filein = pd.read_csv(filename, sep="\t")
        except ValueError as e:
            if "Excel file format cannot be determined" in str(e):
                filein = pd.read_csv(filename, sep="\t")
            else:
                raise
        return filein

    elif filename.endswith(".tsv"):
        return pd.read_csv(filename, sep="\t")
    elif filename.endswith(".csv"):
        return pd.read_csv(filename)
    else:
        raise ValueError("Unknown filetype: %s" % filename)


def freeze(obj):
    """ Turn dicts into frozendict,
        lists into tuples, and sets
        into frozensets, recursively - usefull
        to get a hash value..
    """
    # TODO: combine with ppg.util.freeze
    if hasattr(obj, "__freeze__"):
        return obj.__freeze__()
    else:
        import pypipegraph as ppg

        return ppg.util.freeze(obj)


def parse_a_or_c(ac):
    """parse  an annotator/column combo into a tuple
    anno, column

    Input may be
            a str -> None, column name
            an annotator -> anno, anno.columns[0]
            an (annotator, str) tuple -> anno, str
            an (annotator, int(i)) tuple -> anno, annotator.columns[i]
    """
    from mbf_genomics.annotator import Annotator

    if isinstance(ac, str):
        return (None, ac)
    elif isinstance(ac, Annotator):
        return ac, ac.columns[0]
    elif isinstance(ac, tuple) and len(ac) == 2 and isinstance(ac[0], Annotator):
        if isinstance(ac[1], int):
            return ac[0], ac[0].columns[ac[1]]
        else:
            if not ac[1] in ac[0].columns:
                raise KeyError(
                    "Invalid column name, %s -annotator had %s", (ac[1], ac[0].columns)
                )
            return ac
    elif isinstance(ac, tuple) and len(ac) == 2 and ac[0] is None:
        return ac
    else:
        raise ValueError("parse_a_or_c could not parse %s" % (ac,))


def parse_a_or_c_to_column(k):
    """Parse an annotator + column spec to the column name.
    See parse_a_or_c
    """
    return parse_a_or_c(k)[1]


def parse_a_or_c_to_anno(k):
    """Parse an annotator + column spec to the annotator (or None)
    See parse_a_or_c
    """
    return parse_a_or_c(k)[0]


def parse_a_or_c_to_plot_name(k, default=None):
    """Parse an annotator + column spec to a plot name
    See parse_a_or_c_to_column

    Defaults to column name if no plot_name is defined on annotator

    """
    ac = parse_a_or_c(k)
    if ac[0] is None:
        return k
    return getattr(ac[0], "plot_name", default if default is not None else ac[1])


def find_annos_from_column(k):
    from . import annotator
    import pypipegraph as ppg

    if ppg.inside_ppg():
        if not hasattr(ppg.util.global_pipegraph, "_annotator_singleton_dict"):
            ppg.util.global_pipegraph._annotator_singleton_dict = {}
        singleton_dict = ppg.util.global_pipegraph._annotator_singleton_dict
    else:
        singleton_dict = annotator.annotator_singletons

    res = []
    for anno in singleton_dict["lookup"]:
        if k in anno.columns:
            res.append(anno)
    if res:
        return res
    else:
        raise KeyError("No anno for column '%s' found" % (k,))
