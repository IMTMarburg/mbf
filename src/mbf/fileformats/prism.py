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
