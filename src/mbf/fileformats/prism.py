import pandas as pd
import numpy as np


def read_prism_tables(filename):
    try:
        import BeautifulSoup as bs
        import HTMLParser

        html_decode = HTMLParser.unescape
    except ImportError:
        import bs4 as bs

        # import html
        # bs4 already does this
        html_decode = lambda x: x  # noqa: E731
    with open(filename) as op:
        x = bs.BeautifulSoup(op.read(), "lxml")
    result = []
    for t in x.findAll("table"):
        titles = [html_decode(title.text) for title in t.findAll("title")]
        columns = []
        max_length = 0
        for subcolumn in t.findAll("subcolumn"):
            c = []
            float_count = 0
            for d in subcolumn.findAll("d"):
                dt = html_decode(d.text)
                if dt == "":
                    dt = np.nan
                try:
                    dt = float(dt)
                    float_count += 1
                except ValueError:
                    if dt.count(",") == 1 and dt.count(".") == 0:
                        try:
                            dt = float(dt.replace(",", "."))
                            float_count += 1
                        except ValueError:
                            pass
                c.append(dt)
            if float_count <= 5:
                c = ["" if isinstance(x, float) and np.isnan(x) else x for x in c]
            columns.append(c)
            max_length = max(max_length, len(c))
        for c in columns:
            while len(c) < max_length:
                c.append(np.nan)
        df = pd.DataFrame(
            dict(zip(titles, columns)),
        )[titles]
        result.append(df)
    return result


def read_prism_tables2(filename):
    """When the prism tables are multi leveled?"""
    with open(filename) as op:
        x = bs.BeautifulSoup(op.read(), "lxml")
    result = []
    dfs = []
    for t in x.findAll("table"):
        title = t.findAll("title")[0].text
        x_indices = []
        x_titles = []
        for q in "xcolumn", "xadvancedcolumn":
            for x_column in t.findAll(q):
                title = x_column.find("title").text
                values = [x.text for x in x_column.findAll("d")]
                values = float_or_not(values)
                x_titles.append(title)
                x_indices.append(values)

        columns = []
        column_names = []
        for y in t.findAll("ycolumn"):
            sub_count = len(y.findAll("subcolumn"))
            title = y.find("title").text
            for ii, subcolumn in enumerate(y.findAll("subcolumn")):
                values = [x.text for x in subcolumn.findAll("d")]
                values = float_or_not(values)
                columns.append(values)
                column_names.append((title, str(ii)))
        dd = {k: c for (k, c) in zip(column_names, columns)}
        df = pd.DataFrame(
            dd,
        )
        df.index = pd.MultiIndex.from_tuples(zip(*x_indices), names=x_titles)
        df.name = title
        dfs.append(df)
        break # Todo remove me
        #df.columns = pd.MultiIndex.from_tuples(column_names)
    return dfs
