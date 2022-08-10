"""Genomes from Illuminas Igenomes program"""

from .base import (
    GenomeBase,
    GenomePrebuildMixin,
    class_with_downloads,
    include_in_downloads,
)
from .common import EukaryoticCode
from .filebased import _FileBasedBascis
from mbf_externals.util import (
    download_file_and_gunzip,
    # download_file_and_gzip,
    download_file,
    lazy_property,
    get_page,
    lazy_method,
)
import mbf_externals


@class_with_downloads
class IGenome(GenomePrebuildMixin, _FileBasedBascis):
    def __init__(self, url, archive_date, do_download=True):
        """Get the url from
        https://support.illumina.com/sequencing/sequencing_software/igenome.html
        download it, then have a look for the archive date in <whatever>/Annotation/Archives.
        (the genome build is encoded in the url)

        Here is a list of examples:

        {
            PhiX_RTA: IGenome('http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz', None

        }
        """
        super().__init__()
        self.prebuild_manager = mbf_externals.get_global_manager()

        self.name = self.extract_name_from_url(url)
        self.archive_date = archive_date
        self.url = url
        self.genetic_code = EukaryoticCode
        self.prebuild_prefix = f"igenomes/{self.name}"
        self.download_filename = f"{self.name}.tar.gz"
        self.gtf_filename = 'genes.gtf' if self.archive_date is not None else None
        if do_download:
            self.gene_gtf_dependencies = self.extract_gtf()
            self.cdna_fasta_dependencies = self.create_cdna_from_genome_and_gtf()
            self.download_genome()

    def extract_name_from_url(self, url):
        basename = url[url.rfind("/") + 1 :]
        if not ".tar.gz" in basename:
            raise ValueError(f"could not parse {url}")
        name = basename[: basename.find(".tar.gz")]
        return name

    # todo: merge with ensembl._msg_pack_job

    @include_in_downloads
    @lazy_method
    def download_tar_gz(self):
        fn = self.download_filename
        def do_download(output_dir):
            with open(output_dir / fn, "wb") as op:
                download_file(self.url, op)

        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/download",
            "1",
            [],
            fn,
            do_download,
        )
        self._prebuilds.append(job)  # that was actually missing prior to ppg2.
        return job

    @include_in_downloads
    @lazy_method
    def extract_genome(self):
        fn = "genome.fasta"
        def extract(output_dir):
            import tarfile

            tf = tarfile.open(self.download_tar_gz().find_file(self.download_filename))
            gf = [
                x
                for x in tf.getnames()
                if x.endswith("/Sequence/WholeGenomeFasta/genome.fa")
            ]
            if len(gf) != 1:
                raise ValueError(
                    f"Problem finding WholeGenomeFasta/genome.fasta. Found {gf:r}"
                )
            fh = tf.extractfile(gf[0])
            with open(output_dir / fn, "wb") as op:
                block = fh.read(1024 * 100)
                while block:
                    op.write(block)
                    block = fh.read(1024 * 100)

        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/genome",
            "1",
            [],
            fn,
            extract,
        )
        job.depends_on(self.download_tar_gz())
        self._prebuilds.append(job)  # that was actually missing prior to ppg2.
        return job

    @include_in_downloads
    @lazy_method
    def extract_gtf(self):
        fn = "genes.gtf"
        def extract(output_dir):
            import tarfile

            tf = tarfile.open(self.download_tar_gz().find_file(self.download_filename))
            q = f"/archive-{self.archive_date}/Genes/refGene.txt"
            gf = [
                x
                for x in tf.getnames()
                if x.endswith(q)
            ]
            if len(gf) != 1:
                raise ValueError(
                    f"Problem finding {q}.fasta. Found {gf}. Available {tf.getnames()}"
                )
            fh = tf.extractfile(gf[0])
            with open(output_dir / fn, "wb") as op:
                block = fh.read(1024 * 100)
                while block:
                    op.write(block)
                    block = fh.read(1024 * 100)

        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/genes",
            "1",
            [],
            fn,
            extract,
        )
        job.depends_on(self.download_tar_gz())
        self._prebuilds.append(job)  # that was actually missing prior to ppg2.
        return job

    @include_in_downloads
    def create_cdna_from_genome_and_gtf(self):
        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/cdna",
            "1",
            [],
            "cdna.fasta",
            self._create_cdna_from_genome_and_gtf('cdna.fasta'),
        )
        job.depends_on(
            self.job_transcripts(), self.extract_genome(), self.extract_gtf()
        )
        self._prebuilds.append(job)
        return job

    def create_protein_from_genome_and_gtf(self):
        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/protein",
            "1",
            [],
            "protein.fasta",
            self._create_protein_from_genome_and_gtf('protein.fasta'),
        )
        job.depends_on(self.job_proteins(), self.extract_genome(), self.extract_gtf())
        self._prebuilds.append(job)
        return job
