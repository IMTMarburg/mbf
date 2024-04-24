import scanpy as sc
import shutil
import re
import anndata
import subprocess
from mbf.externals.util import chdir
import os
import pypipegraph2 as ppg
from pathlib import Path
import json


def preprocess_10x_from_star_solo(
    solo_sample, min_genes=200, min_cells=3, dot_size=None, n_genes=6000, pct_mt=15
):
    name = solo_sample.name
    solo_job = solo_sample.load()[0]
    result_dir = Path(f"results/scanpy/preprocessed/{name}")
    result_dir.mkdir(exist_ok=True, parents=True)

    matrix_file = _find_from_job(solo_job, "raw/matrix.mtx.gz")
    barcode_file = _find_from_job(solo_job, "raw/barcodes.tsv.gz")
    features_file = _find_from_job(solo_job, "raw/features.tsv.gz")

    def plot(output_filenames):
        sc.settings.verbosity = (
            3  # verbosity: errors (0), warnings (1), info (2), hints (3)
        )
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor="white")

        # this is straight from the basic tutorial
        # https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

        adata = sc.read_10x_mtx(matrix_file.parent)

        # scanpy plots into cwd...
        os.chdir(result_dir)
        before = adata.shape

        sc.pl.highest_expr_genes(adata, n_top=20, save=True)

        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=3)
        after = adata.shape
        Path("filter_stats.txt").write_text(
            f"""Filtered to
        {after[0]} cells (removed {before[0] - after[0]} cells)
        {after[1]} genes (removed {before[1] - after[1]} genes)
        """
        )

        adata.var["mt"] = adata.var_names.str.startswith(
            "MT-"
        )  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            size=dot_size,
            save="_raw.pdf",
        )

        sc.pl.scatter(
            adata,
            x="total_counts",
            y="pct_counts_mt",
            size=dot_size,
            save="_total_vs_mt.pdf",
        )
        sc.pl.scatter(
            adata,
            x="total_counts",
            y="n_genes_by_counts",
            size=dot_size,
            save="_total_vs_genes.pdf",
        )

        adataf = adata[
            (adata.obs.n_genes_by_counts < n_genes)
            & (adata.obs.pct_counts_mt < pct_mt),
            :,
        ]
        sc.pl.scatter(
            adataf,
            x="total_counts",
            y="pct_counts_mt",
            size=dot_size,
            save="_filtered_vs_mt.pdf",
        )
        sc.pl.scatter(
            adataf,
            x="total_counts",
            y="n_genes_by_counts",
            size=dot_size,
            save="_filtered_vs_genes.pdf",
        )

        adatan = adataf.copy()
        sc.pp.normalize_total(adatan, target_sum=1e4)
        sc.pp.log1p(adatan)
        sc.pp.highly_variable_genes(adatan, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adatan, save=True)

        adatanf = adatan[:, adatan.var.highly_variable]
        sc.pp.regress_out(adatanf, ["total_counts", "pct_counts_mt"])
        sc.pp.scale(adatanf, max_value=10)

        sc.tl.pca(adatanf, svd_solver="arpack")
        sc.pl.pca(adatanf, size=dot_size, save=True)

        from scipy.sparse import csr_matrix

        adatanf.X = csr_matrix(adatanf.X)
        adatanf.write("anndata.h5ad")

    res = ppg.MultiFileGeneratingJob(
        [
            result_dir / x
            for x in [
                "figures/highest_expr_genes.pdf",
                "filter_stats.txt",
                "figures/violin_raw.pdf",
                "figures/scatter_total_vs_genes.pdf",
                "figures/scatter_total_vs_mt.pdf",
                "figures/scatter_filtered_vs_genes.pdf",
                "figures/scatter_filtered_vs_mt.pdf",
                "figures/filter_genes_dispersion.pdf",
                "figures/pca.pdf",
                "anndata.h5ad",
            ]
        ],
        plot,
    ).depends_on(
        matrix_file, barcode_file, features_file
    )  # pull in solo jobs

    res.depends_on_params(
        {
            "min_genes": min_genes,
            "min_cells": min_cells,
            "dot_size": str(dot_size),
            "n_genes": n_genes,
            "pct_mt": "pct_mt",
        }
    )
    res.vid = [solo_sample.vid]
    res.genome = solo_sample.genome
    return res


def combine_and_preprocess_together(
    output_name,
    star_solo_samples,
    map_name=lambda x: x[: x.rfind("_")],
    min_genes=200,
    max_genes=3000,
    min_cells=3,
    pct_mt=20,
    n_top_genes=2000,
    include_raw = False,
):
    """Combine read solo jobs, preprocess them for downstream analysis / cirroculumus

    Result is a hd5ad with normalized expression data,
    and (pca, tsne, umap, leiden, leihoven) on
    (normalized expression data, pearson residuals of 2000 highest variable genes)
    """

    solo_jobs = [solo_sample.load()[0] for solo_sample in star_solo_samples]
    genomes = set((solo_sample.genome for solo_sample in star_solo_samples))
    if len(genomes) != 1:
        raise ValueError(
            "Multiple genomes in combine_and_preprocess_together - not supported"
        )

    result_dir = Path(f"results/scanpy/combined/{output_name}")
    result_dir.mkdir(exist_ok=True, parents=True)

    def combine(output_filenames):
        raw = {}
        for s, j in zip(star_solo_samples, solo_jobs):
            matrix_file = _find_from_job(j, "raw/matrix.mtx.gz")
            path = matrix_file.parent
            adata = sc.read_10x_mtx(path)
            raw[map_name(s.name)] = adata
        if len(raw) != len(star_solo_samples):
            raise ValueError(
                "map_name mapped different samples to same name:", raw.keys()
            )
        combined = anndata.concat(raw, label="dataset")
        combined.obs_names_make_unique()

        adata = combined
        if include_raw:
            adata.write_h5ad(result_dir / f"{output_name}_raw.h5ad")

        before = adata.shape
        # very standard filtering, minimum number of genes per cell,
        # minimum number of cells per gene
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

        # now filter cells that have 'too much' mitochondria
        after = adata.shape
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        # scanpy plots into cwd...
        with chdir(result_dir):
            sc.pl.violin(
                adata,
                ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                jitter=0.4,
                multi_panel=True,
                size=5,
                save="_raw.pdf",
            )

        adataf = adata[
            (adata.obs.n_genes_by_counts < max_genes)
            & (adata.obs.pct_counts_mt < pct_mt),
            :,
        ]
        final = adataf.shape
        Path(result_dir / "filter_stats.txt").write_text(
            f"""Filtered (min_genes / min_cells) to
        {after[0]} cells (removed {before[0] - after[0]} cells)
        {after[1]} genes (removed {before[1] - after[1]} genes)

        Filtered to (mt / max_genes)
        {final[0]} cells (removed {after[0] - final[0]} cells)
        {final[1]} genes (removed {after[1] - final[1]} genes)
        """
        )

        # identify the highly highly variable genes.
        sc.experimental.pp.highly_variable_genes(
            adataf, flavor="pearson_residuals", n_top_genes=n_top_genes
        )

        adatahv = adataf[:, adataf.var["highly_variable"]]
        sc.experimental.pp.normalize_pearson_residuals(adatahv)

        for ad in adataf, adatahv:
            sc.pp.pca(ad, use_highly_variable=True)
            sc.tl.tsne(ad, use_rep="X_pca")
            sc.pp.neighbors(ad)
            sc.tl.umap(ad)
            sc.tl.louvain(ad)
            sc.tl.leiden(ad)

        for k in adatahv.obsm.keys():
            adataf.obsm[k + "_high_variable_genes"] = adatahv.obsm[k]
        for k in ["louvain", "leiden"]:
            adataf.obs[k + "_high_variable_genes"] = adatahv.obs[k]

        sc.pp.log1p(adataf, base=2)

        adataf.write_h5ad(output_filenames[-1])

    output_filenames = [
        result_dir / "figures/violin_raw.pdf",
        result_dir / f"{output_name}.h5ad",
    ]
    if include_raw:
        output_filenames.append(result_dir / f"{output_name}_raw.h5ad")
    res = (
        ppg.MultiFileGeneratingJob(output_filenames, combine)
        .depends_on_params(
            {
                "n_top_genes": n_top_genes,
                "min_genes": min_genes,
                "min_cells": min_cells,
                "pct_mt": pct_mt,
            }
        )
        .self.depends_on(ppg.FunctionInvariant(result_dir / "map_name", map_name))
        .depends_on(solo_jobs)
    )
    res.vid = [x.vid for x in star_solo_samples]
    res.genome = list(genomes)[0]
    return res


class SingleCellForSCB:
    def __init__(self, pre_process_job, name=None, markers={}, uns_files=None):
        """
        @uns is a list of filenames that cirro drops in uns.
        Default is [ "cirro-rank_gene_groups_leiden.json",
            "cirro-rank_gene_groups_louvain.json"].
        Guess it depends on what you did in scanpy
        """
        if uns_files is None:
            uns_files = [
                "cirro-rank_gene_groups_leiden.json",
                "cirro-rank_gene_groups_louvain.json",
            ]

        if isinstance(pre_process_job, (str, Path)):
            h5ad = [Path(pre_process_job)]
            pre_process_job = ppg.global_pipegraph.find_job(h5ad[0])
            assert h5ad[0] in pre_process_job.files
            assert h5ad[0].suffix == ".h5ad"
        else:
            h5ad = [x for x in pre_process_job.files if x.suffix == ".h5ad"]
            if len(h5ad) != 1:
                raise ValueError(
                    "Could not identify (unique) h5ad from that job, found",
                    h5ad,
                    "perhaps pass in filename?",
                )
        self.h5ad = h5ad[0]
        self.prefix = self.h5ad.with_name(self.h5ad.name + "_cirrocumulus")
        self.marker_filename = self.prefix / "markers.json"
        self.output_filenames = [
            self.prefix / (self.prefix.name + ".jsonl"),
            self.prefix / (self.prefix.name + ".jsonl.idx.json"),
        ] + [self.prefix / "uns" / x for x in uns_files]
        if markers:
            self.output_filenames.append(self.marker_filename)
        self.pre_process_job = pre_process_job
        if name is None:
            self.name = self.h5ad.name[: self.h5ad.name.rfind(".")]
        else:
            self.name = name
        self.vid = pre_process_job.vid
        self.genome = pre_process_job.genome
        self.markers = markers

    def load(self):
        def prep(output_filenames):
            if self.prefix.exists():
                shutil.rmtree(self.prefix)
            self.prefix.mkdir()
            if self.markers:
                self.marker_filename.write_text(json.dumps(self.markers))
            cmd = [
                "cirro",
                "prepare_data",
                self.h5ad.absolute(),
                "--format",
                "jsonl",
                "--out",
                str(self.prefix.absolute()),
            ]
            if self.markers:
                cmd.extend(["--markers", str(self.marker_filename.absolute())])
            subprocess.check_call(
                cmd,
                cwd=self.prefix,
                stdout=open(output_filenames[0].with_suffix(".json.stdout"), "wb"),
                stderr=open(output_filenames[0].with_suffix(".json.stderr"), "wb"),
            )
            shutil.move(
                self.prefix.parent / (self.prefix.name + ".jsonl"),
                self.prefix / (self.prefix.name + ".jsonl"),
            )
            shutil.move(
                self.prefix.parent / (self.prefix.name + ".jsonl.idx.json"),
                self.prefix / (self.prefix.name + ".jsonl.idx.json"),
            )

        res = ppg.MultiFileGeneratingJob(self.output_filenames, prep).depends_on(
            self.pre_process_job
        )
        res.depends_on_params({"markers": self.markers})
        return res

    def stats(self):
        adata = anndata.read_h5ad(self.h5ad)
        return {
            "n_cells": adata.shape[0],
            "n_genes": adata.shape[1],
        }


def _find_from_job(job, basename):
    res = []
    for fn in job.files:
        if str(fn).endswith(basename):
            res.append(fn)
    if len(res) == 1:
        return res[0]
    elif len(res) > 1:
        raise ValueError(
            f"found multiple matches for {basename} in {job}. Available: {job.files}"
        )
    elif len(res) == 0:
        raise ValueError(
            f"found no matches for {basename} in {job}. Available: {job.files}"
        )


class CellRangerCount:
    def __init__(
        self, name, sample, transcriptome_directory_path, genome, generate_bam=False
    ):
        """Name is the output name,
        sample a mbf.align.raw.Sample,
        transcriptome_directory_path the ungziped folder from cellranger
        (see https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)

        Todo: the later need to be at least a checksummed, automatic download,
        or better SharedMultiFileGeneratingJob output...

        """
        self.name = name
        ppg.util.assert_uniqueness_of_object(self)
        self.sample = sample
        if not sample.is_paired:
            raise ValueError("Cell ranger needs paired samples")
        self.transcriptome_directory_path = Path(transcriptome_directory_path)
        self.generate_bam = generate_bam
        self.vid = sample.vid[:]
        self.genome = genome

    def load(self):
        result_dir = Path("results/cellranger/")
        result_dir.mkdir(exist_ok=True, parents=True)

        def run(output_filenames):
            input_pairs = self.sample.input_strategy()
            prefixes = set()
            folders = set()
            rex = re.compile(r"^([^_]+)_S\d_L\d{3}_R\d_\d{3}.fastq.gz$")
            for r1, r2 in input_pairs:
                mr1 = rex.findall(r1.name)
                if not mr1:
                    raise ValueError(
                        f"Cell ranger expects filenames like SampleName_S1_L001_R1_001.fastq.gz A, was {r1}"
                    )
                mr2 = rex.findall(r2.name)
                if not mr2:
                    raise ValueError(
                        f"Cell ranger expects filenames like SampleName_S1_L001_R1_001.fastq.gz A, was {r2}"
                    )
                prefixes.add(mr1[0])
                prefixes.add(mr2[0])
                folders.add(str(r1.parent))
                folders.add(str(r2.parent))
            if len(prefixes) > 1:
                raise ValueError(
                    "Multiple fastq prefixes found - not supported (well, technically, if they were the only files in the folder... but you need to extend the detection code for that case"
                )
            if len(folders) > 1:
                raise ValueError("cellranger input must all be in th same folder")
            fastq_folder = Path(list(folders)[0])
            sample = list(prefixes)[0]

            cmd = [
                "cellranger",
                "count",
                f"--id={self.name}",
                f"--transcriptome={self.transcriptome_directory_path.absolute()}",
                f"--fastqs={fastq_folder.absolute()}",
                f"--sample={sample}",
            ]
            if not self.generate_bam:
                cmd.append("--no-bam")
            final_dir = result_dir / self.name
            if final_dir.exists():
                shutil.rmtree(final_dir)
            subprocess.check_call(cmd, cwd=result_dir)

        job = ppg.MultiFileGeneratingJob(
            [
                result_dir / self.name / "outs" / x
                for x in [
                    "molecule_info.h5",
                    "metrics_summary.csv",
                    "filtered_feature_bc_matrix/matrix.mtx.gz",
                    "filtered_feature_bc_matrix/features.tsv.gz",
                    "filtered_feature_bc_matrix/barcodes.tsv.gz",
                    "raw_feature_bc_matrix/matrix.mtx.gz",
                    "raw_feature_bc_matrix/features.tsv.gz",
                    "raw_feature_bc_matrix/barcodes.tsv.gz",
                ]
            ],
            run,
        ).depends_on(self.sample.prepare_input())
        for fn in self.transcriptome_directory_path.glob("**/*"):
            if fn.is_file():
                job.depends_on(ppg.FileInvariant(fn))
        return job


def aggregate_cellranger_jobs(
    output_name,
    cell_ranger_counters,
):
    """Combine cell ranger filtered matrices,

    Raw data is in layer 'raw', there's also sqrt_norm for variance stabilizied analysis downstream
    """
    output_path = Path("results/cellranger/aggr/") / output_name
    output_path.parent.mkdir(exist_ok=True, parents=True)

    def write_csv(output_filename):
        out = "sample_id,molecule_h5,sample\n"
        for cr in cell_ranger_counters:
            h5_file = _find_from_job(cr.load(), "molecule_info.h5")
            out += f"{cr.name},{str(h5_file.absolute())},{cr.name}\n"
        output_filename.write_text(out)

    csv_job = ppg.FileGeneratingJob(
        output_path.parent / f"{output_name}.csv", write_csv
    )

    def run_job(output_filenames):
        if output_path.exists():
            # otherwise, you don't get  the matrix files..
            shutil.rmtree(output_path)
        subprocess.check_call(
            [
                "cellranger",
                "aggr",
                "--id",
                output_name,
                "--csv",
                csv_job.files[0].absolute(),
            ],
            cwd=output_path.parent,
            stdout=open(output_path.parent / f"{output_name}.stdout.txt", "wb"),
            stderr=open(output_path.parent / f"{output_name}.stderr.txt", "wb"),
        )

    res = ppg.MultiFileGeneratingJob(
        [
            output_path.parent / f"{output_name}.stdout.txt",
            output_path.parent / f"{output_name}.stderr.txt",
        ]
        + [
            output_path / "outs" / "count" / x
            for x in [
                # "molecule_info.h5", # not in agg
                # "metrics_summary.csv", # not in agg
                "summary.json",
                "filtered_feature_bc_matrix.h5",
                "filtered_feature_bc_matrix/matrix.mtx.gz",
                "filtered_feature_bc_matrix/features.tsv.gz",
                "filtered_feature_bc_matrix/barcodes.tsv.gz",
                # "raw_feature_bc_matrix/matrix.mtx.gz" # not in agg
                # "raw_feature_bc_matrix/features.tsv.gz", # not in agg
                # "raw_feature_bc_matrix/barcodes.tsv.gz", # not in agg
            ]
        ],
        run_job,
    ).depends_on(csv_job, [x.load() for x in cell_ranger_counters])
    res.vid = [x.vid for x in cell_ranger_counters]
    res.genome = cell_ranger_counters[0].genome

    return res


def cell_ranger_to_scanpy(name, cellranger_job):
    result_dir = Path(f"results/scanpy/cellranger/{name}")
    result_dir.mkdir(exist_ok=True, parents=True)

    h5_file = _find_from_job(cellranger_job, "filtered_feature_bc_matrix.h5")

    def gen(output_filename):
        adata = sc.read_10x_h5(h5_file)

        sc.pp.pca(adata, n_comps=50)

        sc.tl.tsne(adata, use_rep="X_pca")

        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        sc.tl.leiden(adata)

        adata.write_h5ad(output_filename)

    res = ppg.FileGeneratingJob(result_dir / f"{name}.h5ad", gen).depends_on(
        cellranger_job
    )
    res.vid = cellranger_job.vid
    res.genome = cellranger_job.genome
    return res


def _get_scsa_database_path():
    scsa_py_path = subprocess.check_output(["which", "SCSA.py"]).decode().strip()
    return Path(scsa_py_path).resolve().parent.parent / "share" / "whole_v2.db"


def scsa(
    adata,
    is_cancer=False,
    rank_key="rank_genes_groups",
    out_name="scsa",
    fold_change=None,
    p_value=None,
    tissue=None,
    species="Human",
):
    """Perform single cell cell type annotation
    using SCSA ( https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7235421/
    https://github.com/bioinfo-ibms-pumc/SCSA )

    (we have a flake at https://github.com/imTMarburg/flakes/SCSA)

    Adds an obs (out_name) column to the anndate/scanpp adata object, with the best predicted cell type.
    Adds uns (out_name) with the full output, {stdout, stderr, df, cmd}.

    (inplace)!

    use scsa_list_tissues to list available tissues
    """
    import warnings
    import pandas as pd
    import tempfile

    if not species in ("Human", "Mouse"):
        raise ValueError("SCSA only does man and mouse")

    if rank_key not in adata.uns:
        raise ValueError(
            f"adata.uns must have {rank_key}. Perform sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon') first"
        )
    result = adata.uns[rank_key]
    # todo: check if the tissue is valid

    groups = result["names"].dtype.names
    dat = pd.DataFrame(
        {
            group + "_" + key[:1]: result[key][group]
            for group in groups
            for key in ["names", "logfoldchanges", "scores", "pvals"]
        }
    )
    first_gene_name = dat[[x for x in dat.columns if x.endswith("_n")][0]].iloc[0]
    if species == "Human":
        prefix = "ENSG"
    elif species == "Mouse":
        prefix = "ENSMUSG"
    else:
        raise ValueError(
            f"Don't know how to chek that species genes: {species}"
        )  # pragma: no cover
    if not re.match(rf"^{prefix}\d+$", first_gene_name):
        reg = re.compile(rf"{prefix}\d+")
        if reg.search(first_gene_name):

            def extract_ensembl(x):
                try:
                    return reg.findall(x)[0]
                except:  # noqa: E722
                    warnings.warn("No ensembl id found, passing empty to SCSA")
                    return ""

            for c in dat.columns:
                if c.endswith("_n"):
                    dat[c] = [extract_ensembl(x) for x in dat[c]]

        else:
            raise ValueError(
                f"adata.var must contain ensembl genes, {prefix}\\d+.", first_gene_name
            )

    itf = tempfile.NamedTemporaryFile("w", suffix=".csv")
    otf = tempfile.NamedTemporaryFile("w+", suffix=".tsv")
    dat.to_csv(itf.name)
    # nixos specail...
    db = _get_scsa_database_path()

    # fmt: off
    cmd = [
        "SCSA.py",
        "-d", str(db),
        "-i", str(itf.name),
        "-s", "scanpy",
        "-o", str(otf.name),
        "-m", "txt",
        "-g", species,
        "-T", "cancer" if "is_cancer" else "normal",
        # -E # we are using ensembl ids...
        # -f1.5 -p 0.01 maybe later
    ]
    # fmt: on
    if fold_change is not None:
        cmd.extend(["-f", "%.2f" % fold_change])
    if p_value is not None:
        cmd.extend(["-p", "%f" % p_value])
    if tissue:
        cmd.extend(["-k", tissue])
    env = os.environ.copy()
    if "LD_LIBRARY_PATH" in env:  # rpy2 likes to sneak this in, breaking e.g. STAR
        del env["LD_LIBRARY_PATH"]

    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )
    p.check_returncode()

    try:
        otf.seek(0)
        df = pd.read_csv(otf, sep="\t")
    except Exception as e:
        df = str(e)
    uns = {
        "stdout": p.stdout.decode("utf-8"),
        "stderr": p.stderr.decode("utf-8"),
        "df": df,
        "cmd": cmd,
    }
    if isinstance(df, str):
        print("scsa failed. Returning result, no values set on adata")
        return uns
    adata.uns[out_name] = uns
    input_column = adata.uns[rank_key]["params"]["groupby"]
    top = df.groupby("Cluster").head(1).set_index("Cluster")["Cell Type"]
    top.index = pd.Categorical(
        top.index.astype(str), categories=adata.obs[input_column].cat.categories
    )
    adata.obs[out_name] = adata.obs[input_column].replace(top)
    return top


def scsa_list_tissues():
    db = _get_scsa_database_path()
    cmd = [
        "SCSA.py",
        "-d",
        str(db),
        "-i",
        str(db),  # it needs any existing file because bad code...
        "-l",
    ]
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    raw = p.stdout.decode("utf-8")
    result = {}
    species = None
    for line in raw.split("\n"):
        if line.startswith("Species"):
            species = re.findall("Species:([^ ]+)", line)[0]
        elif species and not line.startswith("#") and not line.startswith("-"):
            entries = re.findall(r"\d+:\s+([^0-9]+)", line)
            if not species in result:
                result[species] = []
            result[species].extend([x.strip() for x in entries])
    return result
