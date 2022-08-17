from matplotlib import pyplot as plt
import pypipegraph as ppg
import venn
from functools import partial


def plot_venn(output_prefix, a_dict):
    if hasattr(next(iter(a_dict.values())), "venn_annotator"):
        return plot_venn_from_genes_with_comparisons(output_prefix, a_dict)
    else:
        raise NotImplementedError("Expand!")


def plot_venn_from_genes_with_comparisons(
    output_prefix, a_dict, id_column="gene_stable_id"
):
    if len(a_dict) not in (2, 3):
        raise ValueError("Max support 3 sets currently")

    def plot():
        up = {}
        down = {}
        for name, genes_ddf in sorted(a_dict.items()):
            df = genes_ddf.df
            stable_ids = df[id_column]
            column = genes_ddf.venn_annotator["log2FC"]
            up[name] = set(stable_ids[df[column] > 0])
            down[name] = set(stable_ids[df[column] < 0])
        plt.figure(figsize=(4, 4))
        fixed_venn(up)
        plt.savefig(str(output_prefix) + ".up.png", dpi=72)
        plt.figure(figsize=(4, 4))
        fixed_venn(down)
        plt.savefig(str(output_prefix) + ".down.png", dpi=72)

    return (
        ppg.MultiFileGeneratingJob(
            [str(output_prefix) + ".up.png", str(output_prefix) + ".down.png"], plot
        )
        .depends_on([x.add_annotator(x.venn_annotator) for x in a_dict.values()])
        .depends_on(ppg.ParameterInvariant(output_prefix, id_column))
    )


# the generate_petal_labels in the original venn has a
# bug it won't work with all sets = 0.
# we monkey patch this out here.
def generate_petal_labels(datasets, fmt="{size}"):
    """Generate petal descriptions for venn diagram based on set sizes"""
    datasets = list(datasets)
    n_sets = len(datasets)
    dataset_union = set.union(*datasets)
    universe_size = max(len(dataset_union), 1)
    petal_labels = {}
    for logic in venn._venn.generate_logics(n_sets):
        included_sets = [datasets[i] for i in range(n_sets) if logic[i] == "1"]
        excluded_sets = [datasets[i] for i in range(n_sets) if logic[i] == "0"]
        petal_set = (dataset_union & set.intersection(*included_sets)) - set.union(
            set(), *excluded_sets
        )
        petal_labels[logic] = fmt.format(
            logic=logic,
            size=len(petal_set),
            percentage=(100 * len(petal_set) / universe_size),
        )
    return petal_labels


# it's a bit more involved because they early bind this thing
# but at least we get it cleanly wrapped
# without affecting outside code.
def venn_dispatch(
    data,
    func,
    fmt="{size}",
    hint_hidden=False,
    cmap="viridis",
    alpha=0.4,
    figsize=(8, 8),
    fontsize=13,
    legend_loc="upper right",
    ax=None,
):
    """Check input, generate petal labels, draw venn or pseudovenn diagram"""
    if not venn._venn.is_valid_dataset_dict(data):
        raise TypeError("Only dictionaries of sets are understood")
    if hint_hidden and (func == venn._venn.draw_pseudovenn6) and (fmt != "{size}"):
        error_message = "To use fmt='{}', set hint_hidden=False".format(fmt)
        raise NotImplementedError(error_message)
    n_sets = len(data)
    return func(
        petal_labels=generate_petal_labels(data.values(), fmt),
        dataset_labels=data.keys(),
        hint_hidden=hint_hidden,
        colors=venn._venn.generate_colors(n_colors=n_sets, cmap=cmap, alpha=alpha),
        figsize=figsize,
        fontsize=fontsize,
        legend_loc=legend_loc,
        ax=ax,
    )


fixed_venn = partial(venn_dispatch, func=venn._venn.draw_venn, hint_hidden=False)
