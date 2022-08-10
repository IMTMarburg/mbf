import pypipegraph as ppg
from .base import include_in_downloads, class_with_downloads, DownloadMixin
from mbf_externals.util import lazy_method
import mbf_externals
import pandas as pd

_ensemble_compara_dedup = {}


def EnsemblCompara(version: int, prebuild_manager=None):
    """Singletonize"""
    if prebuild_manager is None:  # pragma: no cover
        prebuild_manager = mbf_externals.get_global_manager()
    if ppg.util.global_pipegraph is not None:
        if not hasattr(ppg.util.global_pipegraph, "_ensemble_compara_dedup"):
            ppg.util.global_pipegraph._ensemble_compara_dedup = {}
        cache = ppg.util.global_pipegraph._ensemble_compara_dedup
    else:
        cache = _ensemble_compara_dedup
    if not version in cache:
        cache[version] = _EnsemblCompara(version, prebuild_manager)
    return cache[version]


@class_with_downloads
class _EnsemblCompara(DownloadMixin):
    def __init__(
        self,
        revision: int,
        prebuild_manager,
        source_species=["Homo_sapiens", "Mus_musculus"],
    ):
        DownloadMixin.__init__(self)
        self.revision = int(revision)
        self.prebuild_manager = prebuild_manager
        self.name = f"Compara_{self.revision}"
        self.source_species = [x.lower() for x in source_species]
        if ppg.inside_ppg():
            ppg.util.assert_uniqueness_of_object(self)

        if hasattr(
            ppg, "is_ppg2"
        ):  # SharedMultiFileGeneartingJobs will sucessfully keep them apart only if necessary
            self.prebuild_prefix = f"ensembl_compara"
        else:
            raise ValueError("not tested on ppg1, sorry")

    @property
    def base_url(self):
        return f"http://ftp.ensembl.org/pub/release-{self.revision}/tsv/ensembl-compara/homologies/"

    @include_in_downloads
    def _download_homology_tables(self):
        jobs = []
        for species in self.source_species:
            for what in "protein", "ncrna":
                jobs.append(
                    self._pb_download_straight(
                        f"{species}/{what}_homologies",
                        species,
                        (
                            rf'href="Compara[.]{self.revision}[.]{what}_default[.]homologies[.]tsv[.]gz"',
                        ),
                        f"compara_{self.revision}_{species}_{what}.tsv.gz",
                        lambda match: "/"
                        + match.strip()[6:-1],  # remove href an quotes
                    )
                )
        return jobs

    def get_homology_table(self, source_species, target_species):
        source_species = source_species.lower()
        target_species = target_species.lower()
        if not source_species in self.source_species:
            raise KeyError(
                f"This compara wasn't setup for {source_species}, add to source_species on EnsemblCompara(...) call"
            )
        result = []
        for what in "protein", "ncrna":
            filename = self.find_file(
                f"compara_{self.revision}_{source_species}_{what}.tsv.gz"
            )
            df = pd.read_csv(filename, sep="\t")
            matching = df[df.homology_species == target_species]
            if len(matching) == 0:
                raise ValueError(
                    f"Target species {target_species} not found in homology table for {self.source_species} - check {filename} for what's in there"
                )
            result.append(matching)
        return pd.concat(result, axis=0)
