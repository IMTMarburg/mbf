import numpy as np
import collections
import pandas as pd
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as pyplot
import skimage

default = object()


def map_to_integers(series, upper, min=None, max=None):
    """Map integers into 0...upper."""
    min = series.min() if min is None else min
    max = series.max() if max is None else max
    zero_to_one = (series - min) / (max - min)
    scaled = zero_to_one * (upper - 1)
    return scaled.astype(int)


def unmap(series, org_series, res):
    """Inverse of map_to_integers"""
    zero_to_one = series / series.max()
    mult = zero_to_one * (org_series.max() - org_series.min())
    shifted = mult + org_series.min()
    return shifted


class ScanpyPlotter:
    """Plotting helpers for anndata objects"""

    def __init__(
        self, ad, cell_type_column="cell_type", colors=default, boundary_resolution=200
    ):
        """
        @ad - ann addata object
        @cell_type_column - which .obs column has your cell type annotation"""
        self.ad = ad
        # wether or not this has "name ENGS" style variable indices
        self.has_name_and_id = ad.var.index.str.contains(" ").any()
        self.cell_type_column = cell_type_column
        if colors is default:
            colors = [
                "#1C86EE",
                "#008B00",
                "#FF7F00",  # orange
                "#4D4D4D",
                "#FFD700",
                "#7EC0EE",
                "#FB9A99",  # lt pink
                "#60D060",  # "#90EE90",
                # "#0000FF",
                "#FDBF6F",  # lt orange
                # "#B3B3B3",
                # "#EEE685",
                "#B03060",
                "#FF83FA",
                # "#FF1493",
                # "#0000FF",
                "#36648B",
                "#00CED1",
                "#00FF00",
                "#8B8B00",
                "#CDCD00",
                "#A52A2A",
            ]
        if isinstance(colors, list):
            cmap = matplotlib.colors.ListedColormap(colors)
        else:
            cmap = colors
        self.cell_type_color_map = cmap

        self.prep_boundaries(boundary_resolution)

    def get_column(self, column):
        """Returns a Series with the data, and the corrected column name"""
        adata = self.ad
        if column in adata.obs:
            pdf = {column: adata.obs[column]}
            column = column
        elif column in adata.var.index:
            pdf = adata[:, adata.var.index == column].to_df()
        else:
            if self.has_name_and_id:
                name_hits = adata.var.index.str.startswith(column + " ")
                if name_hits.sum() == 1:
                    pdf = adata[:, name_hits].to_df()
                    column = pdf.columns[0]
                else:
                    id_hits = adata.var.index.str.endswith(" " + column)
                    if id_hits.sum() == 1:
                        pdf = adata[:, id_hits].to_df()
                        column = pdf.columns[0]
            else:
                raise KeyError("Could not find column %s" % column)
        return pdf[column], column

    def get_column_cell_type(self):
        return self.ad.obs[self.cell_type_column]

    def get_coordinate_dataframe(self, embedding):
        cols = ["x", "y"]
        pdf = (
            pd.DataFrame(self.ad.obsm["X_" + embedding], columns=cols)
            .assign(index=self.ad.obs.index)
            .set_index("index")
        )
        return pdf

    def get_cell_type_categories(self, pdf=None):
        if pdf is None:
            ct = self.get_column_cell_type()
        else:
            ct = pdf["cell_type"]
        if ct.dtype == "category":
            cats = ct.cat.categories
        else:
            cats = sorted(ct.unique())
        return cats

    def prep_boundaries(self, boundary_resolution=200):
        # this image we'll use to find the boundaries
        img = np.zeros((boundary_resolution, boundary_resolution), dtype=np.uint8)
        # and this to determine their colors.
        color_img = np.zeros((boundary_resolution, boundary_resolution), dtype=object)

        pdf = self.get_coordinate_dataframe("umap").assign(
            cell_type=self.get_column_cell_type()
        )

        cats = self.get_cell_type_categories(pdf)

        for cat_no, cat in enumerate(cats):
            sdf = pdf[pdf.cell_type == cat]
            mapped_x = map_to_integers(
                sdf["x"], boundary_resolution, pdf["x"].min(), pdf["x"].max()
            )
            mapped_y = map_to_integers(
                sdf["y"], boundary_resolution, pdf["y"].min(), pdf["y"].max()
            )
            color = self.cell_type_color_map(cat_no)
            for x, y, c in zip(mapped_x, mapped_y, pdf["cell_type"]):
                img[x][y] = True
                color_img[x, y] = color

        flooded = skimage.segmentation.flood(img, (0, 0))
        flooded = skimage.filters.gaussian(flooded)
        bounds = skimage.segmentation.chan_vese(flooded)
        bounds = skimage.segmentation.find_boundaries(bounds)

        # now turn it into something matplotlib can use
        boundary_points = collections.defaultdict(list)

        def search_color(x, y, dist):
            for xi in range(
                max(0, x - dist), min(x + dist + 1, boundary_resolution - 1)
            ):
                for yi in range(
                    max(y - dist, 0), min(y + dist + 1, boundary_resolution - 1)
                ):
                    col = color_img[xi, yi]
                    if col != 0:
                        return col
            return 0

        for x in range(0, boundary_resolution):
            for y in range(0, boundary_resolution):
                if bounds[x][y]:
                    col = color_img[x, y]
                    if col == 0:
                        rdist = 1
                        while rdist < 20 and col == 0:
                            col = search_color(x, y, rdist)
                            rdist += 1
                    if col != 0:
                        boundary_points["x"].append(x)
                        boundary_points["y"].append(y)
                        boundary_points["color"].append(col)
                    else:
                        raise ValueError(
                            "Color was still 0 after looking at 20 cells in each direction? Something is not right"
                        )

        self.boundary_dataframe = pd.DataFrame(
            {
                "x": unmap(
                    pd.Series(boundary_points["x"]), pdf["x"], boundary_resolution
                ),
                "y": unmap(
                    pd.Series(boundary_points["y"]), pdf["y"], boundary_resolution
                ),
                "color": boundary_points["color"],
            }
        )

    def _plot_border_cell_types(
        self, ax, pdf, include_cell_type_legend, border_size, bg_color
    ):
        bdf = self.boundary_dataframe
        ax.scatter(
            bdf["x"],
            bdf["y"],
            color=bdf["color"],
            s=border_size,
            alpha=1,
            edgecolors="none",
            linewidth=0,
            marker="o",
        )

    def _add_legends(
        self,
        fig,
        ax,
        include_cell_type_legend,
        border_cell_types,
        include_color_legend,
        is_numerical,
        cell_type_legend_x_pos,
        subplots_adjust,
    ):
        if (include_cell_type_legend and border_cell_types) or (
            include_color_legend and not is_numerical
        ):
            if cell_type_legend_x_pos is default:
                if include_color_legend:
                    cell_type_legend_x_pos = 1.40
                else:
                    cell_type_legend_x_pos = 1.00

            ax.legend(loc="lower left", bbox_to_anchor=(cell_type_legend_x_pos, 0))

        if subplots_adjust is not default:
            fig.subplots_adjust(right=subplots_adjust)
        else:
            if include_color_legend and include_cell_type_legend:
                fig.subplots_adjust(right=0.65)
                pass
            elif include_cell_type_legend:
                fig.subplots_adjust(right=0.70)

    def plot_cell_density(
        self,
        embedding,
        title=default,
        clip_quantile=0.99,
        border_cell_types=True,
        border_size=15,
        include_color_legend=True,
        include_cell_type_legend=True,
        cell_type_legend_x_pos=default,
        bins=200,
        cmap=default,
        zero_color="#FFFFFF",
        upper_clip_color="#FF0000",
        show_spines=True,
    ):
        pdf = self.get_coordinate_dataframe(embedding).assign(
            cell_type=self.get_column_cell_type()
        )
        hist = np.histogram2d(
            pdf["x"],
            pdf["y"],
            bins=bins,
        )
        vmax_quantile = np.percentile(hist[0], clip_quantile * 100)
        vmax = np.max(hist[0])

        fig, ax = pyplot.subplots(1)
        if border_cell_types:
            self._plot_border_cell_types(
                ax,
                pdf,
                include_cell_type_legend,
                border_size,
                zero_color,
            )
        del hist
        if cmap is default:
            cmap = mcolors.LinearSegmentedColormap.from_list(
                "mine",
                [
                    # "#7f7fFF",
                    "#BFBFFF",
                    "#0000FF",
                ],
                N=256,
            )
            cmap.set_under(zero_color)
            cmap.set_over(upper_clip_color)
        h = ax.hist2d(
            pdf["x"], pdf["y"], bins=bins, vmax=vmax_quantile, cmap=cmap, cmin=1
        )

        def color_map_label(x, pos):
            if x == vmax_quantile:
                if clip_quantile < 1:
                    return ">%.2f" % x
            return "%.2f" % x

        upper = int(np.ceil(vmax))
        step_size = 2
        ticks = list(range(1, upper, step_size))
        if clip_quantile < 1:
            ticks.append(vmax_quantile)
        if include_color_legend:
            cbar = fig.colorbar(
                h[3],
                ax=ax,
                orientation="vertical",
                label="cells per dot",
                extend="both" if clip_quantile < 1 else "min",
                extendrect=True,
                format=matplotlib.ticker.FuncFormatter(color_map_label),
                ticks=ticks,
            )

        self._add_legends(
            fig,
            ax,
            include_cell_type_legend,
            border_cell_types,
            include_color_legend,
            True,
            cell_type_legend_x_pos,
            default,
        )
        # hide the box around the plot
        ax.spines["top"].set_visible(show_spines)
        ax.spines["bottom"].set_visible(show_spines)
        ax.spines["left"].set_visible(show_spines)
        ax.spines["right"].set_visible(show_spines)

        ax.tick_params(
            axis="both",  # changes apply to the x-axis
            which="both",  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False,
        )  # labels along the bottom edge are off

        return fig, ax

    def plot_scatter(
        self,
        embedding,
        gene,
        title=default,
        clip_quantile=0.95,
        border_cell_types=True,
        border_size=15,
        cell_type_legend_x_pos=default,
        plot_zeros=True,
        zero_color="#D0D0D0",
        zero_dot_size=5,
        expression_cmap=default,
        include_color_legend=True,
        include_cell_type_legend=True,
        dot_size=1,
        upper_clip_color="#FF0000",
        subplots_adjust=default,
        plot_data=True,
        bg_color="#FFFFFF",
        plot_categorical_outliers=True,
        categorical_outlier_quantile=0.95,
        anti_overplot=True,
        include_zeros_in_regular_plot=False,
        show_spines=True,
    ):
        expr, expr_name = self.get_column(gene)
        is_numerical = (expr.dtype != "object") and (expr.dtype != "category")

        pdf = (
            self.get_coordinate_dataframe(embedding)
            .assign(expression=expr)
            .assign(cell_type=self.get_column_cell_type())
        )
        fig, ax = pyplot.subplots(layout="tight")
        ax.set_facecolor(bg_color)

        if border_cell_types:
            self._plot_border_cell_types(
                ax,
                pdf,
                include_cell_type_legend,
                border_size,
                bg_color,
            )

        if is_numerical:
            if (
                plot_zeros
            ):  # actually, plot all of them in this color first. That gives you dot sizes to play iwith.
                sdf = pdf  # [pdf["expression"] == 0]
                ax.scatter(
                    sdf["x"],
                    sdf["y"],
                    color=zero_color,
                    s=zero_dot_size,
                    alpha=1,
                    edgecolors="none",
                    linewidth=0,
                    marker=".",
                )
            if include_zeros_in_regular_plot:
                sdf = pdf
            else:  # what you should be doing
                sdf = pdf[pdf["expression"] > 0]
            if anti_overplot:
                sdf = sdf.sort_values("expression")
            expr_min = sdf.expression.min()
            expr_max = sdf.expression.max()
            # add these to the legend
            if expression_cmap is default:
                expression_cmap = mcolors.LinearSegmentedColormap.from_list(
                    "mine",
                    [
                        "#000000",
                        "#0000FF",
                        "#FF00FF",
                    ],
                    N=256,
                )
            cmap_limits = expression_cmap.resampled(256)
            cmap_limits.set_under(zero_color)
            cmap_limits.set_over(upper_clip_color)
            over_threshold = sdf["expression"].quantile(clip_quantile)
            # xdf = sdf[sdf["expression"] >= 0]
            if plot_data:
                plot = ax.scatter(
                    sdf["x"],
                    sdf["y"],
                    c=sdf["expression"],  # .clip(0, upper=over_threshold),
                    cmap=cmap_limits,
                    s=dot_size,
                    alpha=1,
                    vmin=expr_min,
                    vmax=over_threshold,
                )
            else:
                include_color_legend = False

            def color_map_label(x, pos):
                if x == expr_min:
                    return "<%.2f" % x
                elif x == over_threshold:
                    return ">%.2f" % x
                else:
                    return "%.2f" % x

            if include_color_legend:
                ticks = list(range(0, int(np.ceil(expr_max)) + 1))
                if clip_quantile < 1:
                    ticks.append(over_threshold)
                cbar = fig.colorbar(
                    plot,
                    ax=ax,
                    orientation="vertical",
                    label="log2 expression",
                    extend="both",
                    extendrect=True,
                    format=matplotlib.ticker.FuncFormatter(color_map_label),
                    ticks=ticks,
                )
                # this doesn't work
                # cbar.ax.hlines([.110], [0], [1], colors=['red'], linewidth=2)
                if not include_zeros_in_regular_plot:
                    cbar.ax.text(1.50, 0.110, "- 0", ha="center", va="center")

        else:
            if expression_cmap is default:
                cmap = self.cell_type_color_map

            else:
                cmap = expression_cmap
            if pdf["expression"].dtype == "category":
                cats = pdf["expression"].cat.categories
            else:
                cats = pdf["expression"].unique()

            if plot_data:
                for ii, kind in enumerate(cats):
                    sdf = pdf[pdf["expression"] == kind]
                    ax.scatter(
                        sdf["x"],
                        sdf["y"],
                        color=cmap.colors[ii % len(cmap.colors)],
                        s=dot_size,
                        alpha=1,
                        edgecolors="none",
                        linewidth=0,
                        marker=".",
                    )
                    if include_color_legend:
                        ax.scatter(
                            sdf["x"][:0],
                            sdf["y"][:0],
                            color=cmap.colors[ii % len(cmap.colors)],
                            label=kind,
                        )

                # plot the outliers again, so they are on *top* of
                # the regular cell clouds
            if plot_categorical_outliers:
                for ii, kind in enumerate(cats):
                    sdf = pdf[pdf["expression"] == kind]
                    x_center = sdf["x"].mean()
                    y_center = sdf["y"].mean()
                    euclidean_distance = np.sqrt(
                        (sdf["x"] - x_center) ** 2 + (sdf["y"] - y_center) ** 2
                    )
                    threshold = euclidean_distance.quantile(
                        categorical_outlier_quantile
                    )
                    outliers = euclidean_distance > threshold
                    ax.scatter(
                        sdf["x"][outliers],
                        sdf["y"][outliers],
                        color=cmap.colors[ii % len(cmap.colors)],
                        s=dot_size,
                        # color = 'black',
                        alpha=1,
                        edgecolors="none",
                        linewidth=0,
                        marker=".",
                    )
        self._add_legends(
            fig,
            ax,
            include_cell_type_legend,
            border_cell_types,
            include_color_legend,
            is_numerical,
            cell_type_legend_x_pos,
            subplots_adjust,
        )

        ax.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            left=False,
            right=False,
            labelbottom=False,
            labelleft=False,
        )
        # hide x/y axis titles
        ax.set_xlabel("")
        ax.set_ylabel("")

        # hide the box around the plot
        ax.spines["top"].set_visible(show_spines)
        ax.spines["bottom"].set_visible(show_spines)
        ax.spines["left"].set_visible(show_spines)
        ax.spines["right"].set_visible(show_spines)

        ax.tick_params(
            axis="both",  # changes apply to the x-axis
            which="both",  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False,
        )  # labels along the bottom edge are off

        # add a title to the figure
        if title is default:
            title = expr_name
        fig.suptitle(title, fontsize=16)
        return fig, ax


# display(
#     ScanpyPlotter(ad_scrubbed).plot_scatter(
#         "umap",
#         "IL23A",
#         include_cell_type_legend=True,
#         include_color_legend=True,
#         plot_data=True,
#         plot_zeros=True,
#     )
# )
