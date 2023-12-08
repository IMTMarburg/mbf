import numpy as np
import pandas as pd
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as pyplot

default = object()


class ScanpyPlotter:
    """Plotting helpers for anndata objects"""

    def __init__(self, ad, cell_type_column="peng_cell_type"):
        """
        @ad - ann addata object
        @cell_type_column - which .obs column has your cell type annotation"""
        self.ad = ad
        # wether or not this has "name ENGS" style variable indices
        self.has_name_and_id = ad.var.index.str.contains(" ").any()
        self.cell_type_column = cell_type_column

        self.default_colors = [
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

    def get_coordinate_dataframe(self, embedding):
        cols = ["x", "y"]
        pdf = (
            pd.DataFrame(self.ad.obsm["X_" + embedding], columns=cols)
            .assign(index=self.ad.obs.index)
            .set_index("index")
        )
        return pdf

    def plot_scatter(
        self,
        embedding,
        gene,
        title=default,
        clip_quantile=0.95,
        border_celltypes=True,
        border_size=15,
        cell_type_colors=default,  # default | matplotlib.colors.Colormap | List,
        cell_type_legend_x_pos=default,
        plot_zeros=True,
        zero_color="#D0D0D0",
        zero_dot_size=5,
        expression_cmap=default,
        include_color_legend=True,
        include_celltype_legend=True,
        dot_size=1,
        upper_clip_color="#FF0000",
        subplots_adjust=default,
        plot_data=True,
        bg_color="#FFFFFF",
    ):
        expr, expr_name = self.get_column(gene)
        is_numerical = (expr.dtype != "object") and (expr.dtype != "category")

        pdf = (
            self.get_coordinate_dataframe(embedding)
            .assign(expression=expr)
            .assign(cell_type=self.ad.obs[self.cell_type_column])
        )
        fig, ax1 = pyplot.subplots()
        ax1.set_facecolor(bg_color)

        if border_celltypes:
            if cell_type_colors is default:
                cell_type_colors = self.default_colors
            if isinstance(cell_type_colors, list):
                cmap = matplotlib.colors.ListedColormap(cell_type_colors)
            else:
                cmap = cell_type_colors
            if pdf["cell_type"].dtype == "category":
                cats = pdf["cell_type"].cat.categories
            else:
                cats = pdf["cell_type"].unique()

            for ii, cat in enumerate(cats):
                sdf = pdf[pdf["cell_type"] == cat]
                color = cmap.colors[ii % len(cmap.colors)]
                ax1.scatter(sdf["x"], sdf["y"], color=color, s=border_size)
                if include_celltype_legend:
                    ax1.scatter(sdf["x"][:0], sdf["y"][:0], color=color, label=cat)

        if is_numerical:
            if (
                plot_zeros
            ):  # actually, plot all of them in this color first. That gives you dot sizes to play iwith.
                sdf = pdf  # [pdf["expression"] == 0]
                ax1.scatter(
                    sdf["x"],
                    sdf["y"],
                    color=zero_color,
                    s=zero_dot_size,
                    alpha=1,
                    edgecolors="none",
                    linewidth=0,
                    marker=".",
                )
            sdf = pdf[pdf["expression"] > 0].sort_values("expression")
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
                plot = ax1.scatter(
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
                cbar = fig.colorbar(
                    plot,
                    ax=ax1,
                    orientation="vertical",
                    label="log2 expression",
                    extend="both",
                    extendrect=True,
                    format=matplotlib.ticker.FuncFormatter(color_map_label),
                    ticks=[over_threshold] + list(range(0, int(np.ceil(expr_max)) + 1)),
                )
                # this doesn't work
                # cbar.ax.hlines([.110], [0], [1], colors=['red'], linewidth=2)
                cbar.ax.text(1.50, 0.110, "- 0", ha="center", va="center")

        else:
            if expression_cmap is default:
                cmap = matplotlib.colors.ListedColormap(self.default_colors)

            else:
                cmap = expression_cmap
            if plot_data:
                if pdf["expression"].dtype == "category":
                    cats = pdf["expression"].cat.categories
                else:
                    cats = pdf["expression"].unique()
                for ii, kind in enumerate(cats):
                    sdf = pdf[pdf["expression"] == kind]
                    ax1.scatter(
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
                        ax1.scatter(
                            sdf["x"][:0],
                            sdf["y"][:0],
                            color=cmap.colors[ii % len(cmap.colors)],
                            label=kind,
                        )

                # plot the outliers again, so they are on *top* of
                # the regular cell clouds
                for ii, kind in enumerate(pdf["expression"].unique()):
                    sdf = pdf[pdf["expression"] == kind]
                    x_center = sdf["x"].mean()
                    y_center = sdf["y"].mean()
                    euclidean_distance = np.sqrt(
                        (sdf["x"] - x_center) ** 2 + (sdf["y"] - y_center) ** 2
                    )
                    threshold = euclidean_distance.quantile(0.95)
                    outliers = euclidean_distance > threshold
                    ax1.scatter(
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

        if (include_celltype_legend and border_celltypes) or (
            include_color_legend and not is_numerical
        ):
            if cell_type_legend_x_pos is default:
                if include_color_legend:
                    cell_type_legend_x_pos = 1.40
                else:
                    cell_type_legend_x_pos = 1.00

            ax1.legend(loc="lower left", bbox_to_anchor=(cell_type_legend_x_pos, 0))

        if subplots_adjust is not default:
            fig.subplots_adjust(right=subplots_adjust)
        else:
            if include_color_legend and include_celltype_legend:
                fig.subplots_adjust(right=0.65)
                pass
            elif include_celltype_legend:
                fig.subplots_adjust(right=0.70)

        ax1.tick_params(
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
        ax1.set_xlabel("")
        ax1.set_ylabel("")

        # hide the box around the plot
        ax1.spines["top"].set_visible(False)
        ax1.spines["bottom"].set_visible(False)
        ax1.spines["left"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        # add a title to the figure
        if title is default:
            title = expr_name
        fig.suptitle(title, fontsize=16)
        return fig


# display(
#     ScanpyPlotter(ad_scrubbed).plot_scatter(
#         "umap",
#         "IL23A",
#         include_celltype_legend=True,
#         include_color_legend=True,
#         plot_data=True,
#         plot_zeros=True,
#     )
# )
