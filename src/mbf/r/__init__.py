try:
    import rpy2.robjects as ro
    import rpy2.robjects.numpy2ri as numpy2ri  # noqa: F401
    import rpy2.rinterface as rinterface  # noqa: F401
    import rpy2.robjects.pandas2ri as pandas2ri  # noqa: F401
except ImportError:
    import warnings

    warnings.warn("No R available")


__version__ = "0.1"


def convert_dataframe_to_r(obj):
    """Convert a Python DataFRame into int's R equivalent,

    Reimplemented from pandas2ri, but I really don't want to activate the
    automatic.
    """
    with (
        ro.default_converter + numpy2ri.converter + pandas2ri.converter
    ).context() as ctx:
        return ctx.py2rpy(obj)
    # od = {}
    # for name, values in obj.items():
    #     try:
    #         func = pandas2ri.py2rpy.registry[type(values)]
    #         od[name] = func(values)
    #     except Exception as e:  # pragma: no cover - defensive
    #         raise ValueError(
    #             "Error while trying to convert "
    #             'the column "%s". Fall back to string conversion. '
    #             "The error is: %s" % (name, str(e))
    #         )

    # return ro.vectors.DataFrame(od)


def convert_dataframe_from_r(df_r):
    """Take an R dataframe (with colnames and rownames) and turn it into pandas,
    with a reset index."""

    with (
        ro.default_converter + numpy2ri.converter + pandas2ri.converter
    ).context() as ctx:
        return ctx.rpy2py(df_r)
    return pandas2ri.rpy2py_dataframe(df_r)
