import pypipegraph2 as ppg
import re
import pandas as pd
from pathlib import Path
from ..externals import aligners
from ..externals.util import lazy_lookup

barcode_path = Path(__file__).parent / "parsebio_barcode"


def quantify_via_splitpipe(
    sub_libraries, sample_to_barcode_definition, chemistry_version, genome
):
    """
    sub_libraries: a dict: {name: mbf.align.Raw}
    sample_to_barcode_definition: a dictionary mapping sample names to barcodes
        Example:
        {
            "SH01": "A1,A2",
            "SH02": "A3,A4",
            ...
        },
    chemistry_version: a string, 'v2' right now.
    genome: a mbf.genomes.Genome
    """

    sub_jobs = [
        splitpipe.run_on_sublibrary(
            lib,
            sample_to_barcode_definition,
            chemistry_version,
            genome,
        )
        for lib in sub_libraries.values()
    ]

    # todo: Joining the sublibraries.
    # todo: add the 'barcode sequences' from the cell names to the obs. (translate bcid)
    # Do they have different barcodes? Are the same barcodes the same cells?
    # Do I need to differentiate the source?


def quantify_via_starsolo(
    output_name,
    sub_libraries,
    sample_to_barcode_definition,
    chemistry_version,
    genome,
    features="GeneFull",
    strand="Forward",
    umi_dedup="1MM_All",
):
    """
    sub_libraries: a dict of {name: mbf.align.Raw}
    sample_to_barcode_definition: a dictionary mapping sample names to barcodes
        Example:
        {
            "SH01": "A1,A2",
            "SH02": "A3,A4",
            ...
        },
    chemistry_version: a string, 'v2' right now.
    genome: a mbf.genomes.Genome

    Returns jobs {'alignment': {}, 'h5ad': job}
    """
    import anndata
    assert chemistry_version == "v2"
    if chemistry_version == "v2":
        barcode_list = "n24_v4", "v1", "v1"
    else:
        raise ValueError(
            "Unknown/undefined chemistry version, update code to use right barcode list",
            chemistry_version,
        )

    assert features in "GeneFull", "Gene"
    star_solo = aligners.STARSolo()
    for lane in sub_libraries.values():
        if not lane.is_paired:
            raise ValueError("Unpaired lane?, not parsebio", lane)

    alignments = {
        k: v.align(
            star_solo,
            genome,
            {
                "soloType": "CB_UMI_Complex",
                "soloCBposition": [
                    "0_78_0_85",
                    "0_48_0_55",
                    "0_10_0_17",
                ],  # star seems right inclusive in it's coordinates.
                "soloUMIposition": "0_0_0_10",
                "soloCBmatchWLtype": "EditDist_2",
                "soloCellFilter": "None",  # do not filter for now, we still need to combine by barcode.
                "cell_barcode_whitelist": get_starsolo_whitelist_jobs(barcode_list),
                "soloFeatures": features,
                "soloStrand": strand,
                "soloUMIdedup": umi_dedup,
            },
            name=v.name + f"_{features}_{strand}_{umi_dedup}",
        )
        for k, v in sub_libraries.items()
    }
    alignment_jobs = {k: v.load() for (k, v) in alignments.items()}
    matrix_files = {
            k: [x for x in 
                v[0] # that's a alignment, bai generation tuple...
                .filenames if x.name == "matrix.mtx.gz"][0]
            for (k,v) in alignment_jobs.items()}
    feature_files = {
            k: [x for x in 
                v[0] # that's a alignment, bai generation tuple...
                .filenames if x.name == "features.tsv.gz"][0]
            for (k,v) in alignment_jobs.items()}
    barcode_files = {
            k: [x for x in 
                v[0] # that's a alignment, bai generation tuple...
                .filenames if x.name == "barcodes.tsv.gz"][0]
            for (k,v) in alignment_jobs.items()}
    assert matrix_files
    assert len(matrix_files) == len(feature_files) == len(barcode_files)


    well_to_sample = interpret_sample_to_barcode_definition(
        sample_to_barcode_definition
    )

    def annotate_and_combine(output_filename):
        replacements = load_barcode_replacements(barcode_list[0])

        ads = []
        for k in sub_libraries.keys():
            ad = load_star_data(
                matrix_files[k].parent 
            )
            rows_to_add_to, rows_to_add_from = determine_replacements(ad, replacements)
            # now add and generate a new ad
            out = ad.X[rows_to_add_to, :]
            add = ad.X[rows_to_add_from, :]
            out += add
            obs = ad.obs.iloc[rows_to_add_to]
            ad2 = anndata.AnnData(out, obs, ad.var)
            ad2.obs["sample"] = [well_to_sample[x[:x.find("_")] for x in ad2.obs.index]
            ads.append(ad2)
        ad = anndata.concat(
            ads, merge="same", label="sublibrary", keys=list(sub_libraries.keys())
        )
        ad.write_h5ad(output_filename)

    combo_job = ppg.FileGeneratingJob(
        Path("results/StarSolo") / output_name / "output_name.h5ad",
        annotate_and_combine,
    )
    combo_job.depends_on(matrix_files.values())
    combo_job.depends_on(feature_files.values())
    combo_job.depends_on(barcode_files.values())
    combo_job.depends_on_file(barcode_path / f"bc_data_{barcode_list[0]}.csv")
    combo_job.depends_on_func(load_star_data)
    combo_job.depends_on_func(load_barcode_replacements)
    combo_job.depends_on_func(determine_replacements)
    combo_job.depends_on_params(
        (
            well_to_sample,
            tuple(sorted(sub_libraries.keys())),
            tuple(sorted([x.name for x in sub_libraries.values()])),
        ),
    )

    return {"alignment": alignment_jobs}


def determine_replacements(ad, replacements):
    rows_to_add_to = []
    rows_to_add_from = []
    for ii, barcode in enumerate(ad.obs.index):
        offset = barcode.find("_")
        bc1 = barcode[:offset]
        if bc1 in replacements:
            partner = replacements[bc1] + barcode[offset:]
            try:
                partner_ii = ad.obs.index.get_loc(partner)
                rows_to_add_to.append(ii)
                rows_to_add_from.append(partner_ii)
            except KeyError:
                raise ValueError(
                    "Replacement barcode not found. I mean it could happen with very low sequenced cells? Not anticipating it for now."
                )
    return (rows_to_add_to, rows_to_add_from)


def load_barcode_replacements(barcode_name):
    """Load the barcode(1) (R) -> barcode(1) R replacement list"""
    barcode_info = pd.read_csv(barcode_path / f"bc_data_{barcode_name}.csv", sep=",")
    # We combine T and R into T
    # to be identical to splitpipe's approach.
    pairs_by_well = {}
    for _, row in barcode_info.iterrows():
        well = row["well"]
        if not well in pairs_by_well:
            pairs_by_well[well] = [None, None]
        if row["stype"] == "T":
            pairs_by_well[well][0] = row["sequence"]
        else:
            pairs_by_well[well][1] = row["sequence"]
    replacements = {x[1]: x[0] for x in pairs_by_well.values()}
    assert len(replacements) == len(barcode_info) / 2
    return replacements


def interpret_sample_to_barcode_definition(sample_to_barcode_definition):
    """Convert barcode to sample definitions like splitpipe.

    For now, we only take straight barcodes, not the block/range definition
    splitpipe does.

    """
    well_to_sample = {}
    for sample, barcodes in sample_to_barcode_definition.items():
        for barcode in barcodes.split(","):
            if not re.match("[A-H][0-9]{1,2}", barcode):
                raise ValueError(
                    "Invalid barcode - not taking blocks/ranges yet", barcode
                )
            well_to_sample[barcode] = sample
    return well_to_sample


def load_star_data(input_path):
    import anndata

    input_path = Path(input_path)
    ad = anndata.read_mtx(input_path / "matrix.mtx.gz")
    ad = ad.transpose()
    features = pd.read_csv(input_path / "features.tsv.gz", sep="\t", header=None)
    barcodes = pd.read_csv(input_path / "barcodes.tsv.gz", sep="\t", header=None)[0]
    ad.obs_names = barcodes
    ad.var_names = features[1] + " " + features[0]
    return ad


def convert_to_starsolo_whitelist(input_filename):
    input_filename = Path(input_filename)
    output_filename = (
        Path("cache") / "starsolo_barcode_whitelists" / input_filename.name
    )
    output_filename.parent.mkdir(exist_ok=True, parents=True)

    def do_conversion(of):
        input = input_filename.read_text().strip().split("\n")[1:]
        output = [x.split(",")[1] for x in input]
        of.write_text("\n".join(output))

    return (
        ppg.FileGeneratingJob(output_filename, do_conversion)
        .depends_on_file(input_filename)
        .self
    )


@lazy_lookup
def get_starsolo_whitelist_jobs(barcode_names):
    return [
        convert_to_starsolo_whitelist(barcode_path / f"bc_data_{barcode_name}.csv")
        for barcode_name in barcode_names
    ]
