import numpy as np
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
            # "#E31A1C",  # red
            "#008B00",
            # "#6A3D9A",  # purple
            "#FF7F00",  # orange
            "#4D4D4D",
            "#FFD700",
            "#7EC0EE",
            "#FB9A99",  # lt pink
            "#90EE90",
            # "#0000FF",
            "#FDBF6F",  # lt orange
            # "#B3B3B3",
            # "#EEE685",
            "#B03060",
            "#FF83FA",
            "#FF1493",
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
            pdf = [adata.obs[column]]
            column = 0
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
        cell_type_colors=default,
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
    ):
        expr, expr_name = self.get_column(gene)
        pdf = (
            self.get_coordinate_dataframe(embedding)
            .assign(expression=expr)
            .assign(cell_type=self.ad.obs[self.cell_type_column])
        )
        fig, ax1 = pyplot.subplots()

        if border_celltypes:
            if cell_type_colors is default:
                cell_type_colors = self.default_colors
            cell_type_cmap = {
                cat: cell_type_colors[ii]
                for ii, cat in enumerate(sorted(pdf["cell_type"].unique()))
            }
            for cat, color in cell_type_cmap.items():
                sdf = pdf[pdf["cell_type"] == cat]
                ax1.scatter(sdf["x"], sdf["y"], color=color, s=border_size)
                if include_celltype_legend:
                    ax1.scatter(sdf["x"][:0], sdf["y"][:0], color=color, label=cat)

        if plot_zeros:
            sdf = pdf[pdf["expression"] == 0]
            ax1.scatter(sdf["x"], sdf["y"], color=zero_color, s=zero_dot_size)

        sdf = pdf[pdf["expression"] > 0]
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
                # sdf_range["expression"].max(),
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
            fig.colorbar(
                plot,
                ax=ax1,
                orientation="vertical",
                label="log2 expression",
                extend="both",
                extendrect=True,
                format=matplotlib.ticker.FuncFormatter(color_map_label),
                ticks=[0, sdf_range.expression.min(), over_threshold]
                + list(range(0, int(np.ceil(sdf_range.expression.max())) + 1)),
            )

        if include_celltype_legend:
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
