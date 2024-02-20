from pathlib import Path
from abc import ABC, abstractmethod
import pandas as pd
from dppd import dppd
import pysam
from .common import reverse_complement, df_to_rows
from .gene import Gene, Transcript
from mbf.externals.util import (
    lazy_method,
    download_file_and_gunzip,
    get_page,
    download_file,
)
import weakref
import numpy as np
import pypipegraph as ppg
import re

dp, X = dppd()


def include_in_downloads(func):
    """A decorator to collect the download funcs"""
    func._include_in_downloads = True
    return func


def class_with_downloads(cls):
    """class decoractor
    to collect functions tagget with @include_in_downloads.

    Use with DownloadMixin
    """
    cls._download_methods = []
    for f in cls.__dict__.items():
        if hasattr(f[1], "_include_in_downloads"):
            cls._download_methods.append(f[1])
    return cls


def ReadOnlyPropertyWithFunctionAccess(func):
    """With normal property, you can not (easily) retrieve
    the function. This will return the value of the func
    if you do x.prop and the func itself if you do type(x).prop
    """

    class Property:
        def __get__(self, inst, instcls):
            if inst is None:
                # instance attribute accessed on class, return self
                return func
            else:
                return func(inst)

    return Property()


class MsgPackProperty:
    """
    a message pack property is a property x_y that get's
    calculated by a method _prepare_x_y
    and automatically stored/loaded by a caching job
    as a file. (The file used to be msgpack, nowadays it's parquet, msgpack
                is no longer supported by pandas, and the mbf-msg-pack has been
                bitrotting ever since)
    the actual job used depends on the GenomeBase subclass

    The dependency_callback get's called with the GenomeBase subclass
    instance and can return dependencys for the generated job

    The object has three members afterwards:
        x_y -> get the value returned by _prepare_x_y (lazy load)
        _prepare_x_y -> that's the one you need to implement,
                    it's docstring is copied to this propery
        job_y -> the job that caches _prepare_x_y() results
    Optionally, impement
        _fix_after_load_x_y -> this is what

    """

    def __init__(self, dependency_callback=None, files_to_invariant_on_callback=None):
        self.dependency_callback = dependency_callback
        self.files_to_invariant_on_callback = files_to_invariant_on_callback


def msgpack_unpacking_class(cls):
    msg_pack_properties = []
    for d in list(cls.__dict__):
        v = cls.__dict__[d]
        if isinstance(v, MsgPackProperty):
            if not "_" in d:
                raise NotImplementedError(
                    "Do not know how to create job name for msg_pack_properties that do not containt  _"
                )
            msg_pack_properties.append(d)
            job_name = "job_" + d[d.find("_") + 1 :]
            filename = d + ".parquet"
            calc_func = getattr(cls, f"_prepare_{d}")

            def load(self, d=d, filename=filename, job_name=job_name):
                if not hasattr(self, "_" + d):
                    fn = self.find_file(filename)
                    if not fn.exists():
                        raise ValueError(
                            f"{d} accessed before the respecting {job_name} call - {fn.absolute()} did not exist"
                        )
                    df = pd.read_parquet(fn)
                    if hasattr(self, "_fix_after_load_" + d):
                        df = getattr(self, "_fix_after_load_" + d)(df)
                    setattr(self, "_" + d, df)
                return getattr(self, "_" + d)

            p = property(load)
            p.__doc__ == calc_func.__doc__
            setattr(cls, d, p)
            if not hasattr(cls, job_name):

                def gen_job(
                    self,
                    d=d,
                    filename=filename,
                    calc_func=calc_func,
                    dependency_callback=v.dependency_callback,
                    files_to_invariant_on_callback=v.files_to_invariant_on_callback,
                ):
                    if files_to_invariant_on_callback:
                        files_to_invariant_on = files_to_invariant_on_callback(self)
                    else:
                        files_to_invariant_on = []
                    j = self._msg_pack_job(
                        d, filename, calc_func, files_to_invariant_on
                    )
                    if j is not None:
                        j.depends_on(dependency_callback(self))
                    return j

                setattr(cls, job_name, gen_job)
            else:  # pragma: no cover
                pass
    if hasattr(cls, "_msg_pack_properties"):
        msg_pack_properties.extend(cls._msg_pack_properties)
    cls._msg_pack_properties = msg_pack_properties
    return cls


class DownloadMixin:
    def __init__(self):
        self._prebuilds = []
        self._download_jobs = []
        self._filename_lookups = []

    def find_file(self, name):
        if name in self._filename_lookups:
            return self._filename_lookups[name]
        for job in self._prebuilds:
            if hasattr(job, "find_file"):  # and not hasattr(ppg, 'is_ppg2'):
                try:
                    return job.find_file(name)
                except KeyError:
                    pass
            else:
                for j in job:
                    for f in j.files:
                        if Path(f).name == name:
                            return Path(f)
        # now search for undeclared, but created files
        # mostly for aligners, where we only track the sentinels, not the index
        # files
        # ppg2: if aligners don't cleare their files, blame them
        if not hasattr(ppg, "is_ppg2"):
            for job in self._prebuilds:
                if hasattr(job, "name_file"):  # pragma: no branch
                    if job.name_file(name).exists():
                        return job.name_file(name)
        raise OSError(f"File not found: {name}")

    def find_prebuild(self, name):
        """Find which prebuild created the file named @name.
        Must be in the list of job.files"""

        for job in self._prebuilds:
            if hasattr(job, "find_file") and not hasattr(ppg, "is_ppg2"):
                # ppg2 you have no target_file at the point you 're trying to find a prebuild...
                try:
                    job.find_file(name)
                    return job
                except KeyError:
                    pass
            else:
                for j in job:
                    for f in j.files:
                        if Path(f).name == name:
                            return job
        raise OSError(f"File not found: {name}")

    @lazy_method
    def download(self):
        """All the jobs needed to download the genome and prepare it for usage"""
        result = []
        for method in self.__class__._download_methods:
            j = method(self)
            if isinstance(j, list):
                if j is not None:  # pragma: no branch
                    result.extend(j)
            elif j is not None:  # pragma: no branch
                result.append(j)
            # for j in result:
            # if isinstance(j, list):
            # raise ValueError(method)
        for j in self._download_jobs:
            if j is not None:  # pragma: no branch
                result.append(j)
        for j in result:
            if not j in self._prebuilds:  # pragma: no branch
                self._prebuilds.append(j)
        if hasattr(self.__class__, "_msg_pack_properties"):
            for msg_pack_prop in self.__class__._msg_pack_properties:
                job_name = "job" + msg_pack_prop[msg_pack_prop.find("_") :]
                j = getattr(self, job_name)()
                self._prebuilds.append(j)
        return result

    def _pb_download(
        self,
        pb_name,
        url,
        regexps,
        output_filename,
        download_func,
        match_transformer=lambda x: x,
        additional_input_files=[],
    ):
        """regexps may be multiple - then the first one matching is used"""

        def do_download(output_path):
            real_url = self.base_url + url
            raw = get_page(real_url)
            if not raw:  # pragma: no cover
                raise ValueError("Retrieving url failed: %s" % real_url)
            for aregexps in regexps:
                matches = re.findall(aregexps, raw)
                if len(matches) == 1:
                    Path(str(output_path / output_filename) + ".url").write_text(
                        (real_url + match_transformer(matches[0]))
                    )
                    download_func(
                        real_url + match_transformer(matches[0]),
                        output_path / output_filename,
                    )
                    break
                else:
                    print(aregexps, len(matches), "found")
            else:
                import sys

                sys.stderr.write(raw)
                raise ValueError(  # pragma: no cover - defensive
                    "Found either too few or too many for every regexps. \nRaw was (see stderr)"
                )

            if Path(output_filename).suffix == ".fasta":
                import pysam

                pysam.faidx(str((output_path / output_filename).absolute()))

        output_filenames = [output_filename]
        for k in output_filenames[:]:
            if Path(k).suffix == ".fasta":
                output_filenames.append(Path(k).with_suffix(".fasta.fai"))

        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/{pb_name}",
            "1",
            [] + additional_input_files,
            output_filenames,
            do_download,
        )
        if hasattr(self, "_pb_find_server"):
            job.depends_on(self._pb_find_server())
        self._prebuilds.append(job)  # that was actually missing prior to ppg2.
        return job

    def _pb_download_straight(
        self,
        pb_name,
        url,
        regexps,
        output_filename,
        match_transformer=lambda x: x,  # pragma: no cover
    ):
        def df(url, filename):
            with open(filename, "wb") as op:
                download_file(url, op)

        # create the jobs only once, or you run into problems with ppg2s more sensitive
        # closure omparisobic
        if not hasattr(self, "_pb_download_straight_jobs"):
            self._pb_download_straight_jobs = {}
        if not output_filename in self._pb_download_straight_jobs:
            self._pb_download_straight_jobs[output_filename] = self._pb_download(
                pb_name, url, regexps, output_filename, df, match_transformer
            )
        return self._pb_download_straight_jobs[output_filename]

    def _pb_download_and_gunzip(self, pb_name, url, regexps, output_filename):
        return self._pb_download(
            pb_name, url, regexps, output_filename, download_file_and_gunzip
        )


@msgpack_unpacking_class
class GenomeBase(ABC, DownloadMixin):
    def __init__(self):
        DownloadMixin.__init__(self)

    @abstractmethod
    def _msg_pack_job(
        self, property_name, filename, callback_function, files_to_invariant_on
    ):
        raise NotImplementedError  # pragma: no cover

    def download_genome(self):
        return self.download()

    @lazy_method
    def get_chromosome_lengths(self):
        """Return a dict name -> length for the primary assembly"""
        f = pysam.FastaFile(str(self.find_file("genome.fasta")))
        return dict(zip(f.references, f.lengths))

    @lazy_method
    def get_true_chromosomes(self):
        """Get the names of 'true' chromosomes, ie. no scaffolds/contigs
        in genomes that have chromosomes, otherwise all"""
        return list(self.get_chromosome_lengths().keys())

    def get_genome_sequence(self, chr, start, stop):
        f = pysam.FastaFile(str(self.find_file("genome.fasta")))
        return f.fetch(chr, start, stop)

    def get_cdna_sequence(self, transcript_stable_id):
        with pysam.FastaFile(str(self.find_file("cdna.fasta"))) as f:
            return f.fetch(transcript_stable_id)

    def get_cds_sequence(self, protein_id, protein_info=None):
        """Get the coding sequence (rna) of a protein"""
        if protein_info is None:
            protein_info = self.df_proteins.loc[protein_id]
        elif protein_info.name != protein_id:
            raise ValueError("protein_id != protein_info['protein_id']")
        cdna = ""
        chr = protein_info["chr"]
        for start, stop in protein_info["cds"]:
            cdna += self.get_genome_sequence(chr, start, stop)
        if protein_info["strand"] not in (1, -1):  # pragma: no cover
            raise ValueError(f'{protein_info["strand"]} was not 1/-1')
        if protein_info["strand"] == -1:
            cdna = reverse_complement(cdna)
        return cdna

    def get_protein_sequence(self, protein_id):
        """Get the AA sequence of a protein"""
        with pysam.FastaFile(str(self.find_file("pep.fasta"))) as f:
            return f.fetch(protein_id)

    def get_additional_gene_gtfs(self):
        return []

    def get_gtf(self, features=[]):
        import mbf_gtf

        filenames = [self.find_file("genes.gtf")]
        filenames.extend(self.get_additional_gene_gtfs())
        dfs = {}
        for gtf_filename in filenames:
            if gtf_filename is None:
                pass
            else:
                r = mbf_gtf.parse_ensembl_gtf(str(gtf_filename), list(features))
                for k, df in r.items():
                    if not k in dfs:
                        dfs[k] = []
                    dfs[k].append(df)
        for k in features:
            if not k in dfs:
                dfs[k] = [pd.DataFrame({})]
        result = {k: pd.concat(dfs[k], sort=False) for k in dfs}
        return result

    @property
    def genes(self):
        """a  dictionary gene_stable_id -> gene.Gene"""
        if not hasattr(self, "_genes"):
            self.build_genes_and_transcripts()
        return self._genes

    @property
    def transcripts(self):
        """a  dictionary transcript_stable_id -> gene.Transcript"""
        if not hasattr(self, "_transcripts"):
            self.build_genes_and_transcripts()
        return self._transcripts

    def name_to_gene_ids(self, name):
        if not hasattr(self, "_name_to_gene_lookup"):
            lookup = {}
            for a_name, stable_id in zip(self.df_genes["name"], self.df_genes.index):
                a_name = a_name.upper()
                if not a_name in lookup:
                    lookup[a_name] = set([stable_id])
                else:
                    lookup[a_name].add(stable_id)
            self._name_to_gene_lookup = lookup
        try:
            return set(self._name_to_gene_lookup[name.upper()])
        except KeyError:
            return set()
        # return set(self.df_genes.index[self.df_genes.name.str.upper() == name.upper()])

    def build_genes_and_transcripts(self):
        genes = {}
        for tup in self.df_genes.itertuples():
            g = Gene(
                tup[0],
                tup.name,
                tup.chr,
                tup.start,
                tup.stop,
                tup.strand,
                tup.biotype,
                transcripts=[],
                genome=weakref.proxy(self),
            )
            genes[tup[0]] = g
        transcripts = {}
        for tup in self.df_transcripts.itertuples():
            g = genes[tup.gene_stable_id]
            t = Transcript(
                tup[0],
                tup.gene_stable_id,
                tup.name,
                tup.chr,
                tup.start,
                tup.stop,
                tup.strand,
                tup.biotype,
                tuple(tup.exons),
                tup.exon_stable_ids,
                weakref.proxy(g),
                genome=weakref.proxy(self),
            )
            transcripts[tup[0]] = t
            g.transcripts.append(t)
        self._genes = genes
        self._transcripts = transcripts

    @ReadOnlyPropertyWithFunctionAccess
    def df_exons(self):
        """a dataframe of all exons (on canonical chromosomes - ie those in get_chromosome_lengths())"""
        res = {
            "chr": [],
            "start": [],
            "stop": [],
            "transcript_stable_id": [],
            "gene_stable_id": [],
            "strand": [],
        }
        canonical_chromosomes = self.get_chromosome_lengths()
        for tr in self.transcripts.values():
            if not tr.chr in canonical_chromosomes:  # pragma: no cover
                continue
            for start, stop in tr.exons:
                res["chr"].append(tr.chr)
                res["start"].append(start)
                res["stop"].append(stop)
                res["transcript_stable_id"].append(tr.transcript_stable_id)
                res["gene_stable_id"].append(tr.gene_stable_id)
                res["strand"].append(tr.strand)
        return pd.DataFrame(res)

    def _prepare_df_genes(self):
        """Return a DataFrame with  gene information:
        gene_stable_id
        name
        chr
        start
        stop
        strand
        tss
        tes
        biotype
        """
        gtf = self.get_gtf(["gene", "transcript"])
        genes = gtf["gene"]
        transcripts = gtf["transcript"]
        if len(genes) == 0:  # a genome without gene information
            return pd.DataFrame(
                {
                    "gene_stable_id": [],
                    "name": [],
                    "chr": [],
                    "start": [],
                    "stop": [],
                    "strand": [],
                    "tss": [],
                    "tes": [],
                    "biotype": [],
                }
            )
        elif len(transcripts) == 0:  # pragma: no cover
            raise ValueError(
                "Genome with gene but no transcript information "
                "not supported: len(genes) %i, len(transcripts) %i"
                % (len(genes), len(transcripts))
            )

        transcripts = transcripts.set_index("gene_id").sort_values(
            ["seqname", "start", "end"]
        )
        genes = (
            dp(genes)
            .transassign(
                gene_stable_id=X.gene_id,
                name=list(
                    X.gene_name if hasattr(X, "gene_name") else X.gene_id
                ),  # this makes sure we have a str(object) column in the dataframe
                # which triggers msgpack not to mess up our tuple columns.
                chr=pd.Categorical(X.seqname),
                start=X.start,
                stop=X.end,
                strand=X.strand,
                tss=(genes.start).where(
                    genes.strand == 1, genes.end
                ),  # note that this is pd.Series.where, which is (cond, other)
                tes=(genes.end).where(genes.strand == 1, genes.start),
                biotype=pd.Categorical(X.gene_biotype),
            )
            .sort_values(["chr", "start"])
            .set_index("gene_stable_id")
            .pd
        )
        if not genes.index.is_unique:
            raise ValueError("gene_stable_ids were not unique")
        tr = {}
        for gene_stable_id, transcript_stable_id in transcripts[
            "transcript_id"
        ].items():
            if not gene_stable_id in tr:
                tr[gene_stable_id] = []
            tr[gene_stable_id].append(transcript_stable_id)
        genes = genes.assign(
            transcript_stable_ids=pd.Series(list(tr.values()), index=list(tr.keys()))
        )
        self.sanity_check_genes(genes)
        return genes

    def sanity_check_genes(self, df_genes):
        strand_values = set(df_genes.strand.unique())
        if strand_values.difference([1, -1]):  # pragma: no cover
            # this is currently already being handled by the gtf parser - defensive
            raise ValueError(f"Gene strand was outside of 1, -1: {strand_values}")
        wrong_order = df_genes["start"] > df_genes["stop"]
        if wrong_order.any():
            raise ValueError("start > stop %s" % df_genes[wrong_order].head())

    def _prepare_df_transcripts(self):
        """Get a DataFrame with all the transcript information
        transcript_stable_id (index),
        gene_stable_id,
        name,
        chr, start, stop, strand,
        tss, tes
        biotype,
        exons  - list (start, stop)
        """
        gtf = self.get_gtf(["transcript", "exon"])
        transcripts = gtf["transcript"]
        exons = gtf["exon"]
        if len(transcripts) == 0:
            df = (
                pd.DataFrame(
                    {
                        "transcript_stable_id": [],
                        "gene_stable_id": [],
                        "name": [],
                        "chr": [],
                        "start": [],
                        "stop": [],
                        "strand": [],
                        "biotype": [],
                        "exons": [],
                        "exon_stable_ids": [],
                        "translation_start": [],
                        "translation_start_exon": [],
                        "translation_stop": [],
                        "translation_stop_exon": [],
                        "protein_id": [],
                    }
                )
                .set_index("transcript_stable_id")
                .sort_values(["chr", "start"])
            )
            return df
        all_exons = exons.set_index("transcript_id").sort_values("start")

        result = (
            dp(transcripts)
            .transassign(
                transcript_stable_id=X.transcript_id,
                gene_stable_id=X.gene_id,
                name=X.transcript_name
                if hasattr(X, "transcript_name")
                else X.transcript_id,
                chr=pd.Categorical(X.seqname),
                start=X.start,
                stop=X.end,
                tss=(X.start).where(
                    X.strand == 1, X.end
                ),  # note that this is pd.Series.where, which is (cond, other)
                tes=(X.end).where(X.strand == 1, X.start),
                strand=X.strand,
                biotype=pd.Categorical(X.transcript_biotype),
                msg_pack_fix=""  # stupid msg_pack writer will mess up the exon tuples
                # if it has no str-object columns in the datafram.
            )
            .set_index("transcript_stable_id")
            .pd
        )

        if not result.index.is_unique:
            raise ValueError("transcript_stable_ids were not unique")
        result_exons = {}
        result_exon_ids = {}
        for transcript_stable_id, estart, estop, eid in zip(
            all_exons.index, all_exons["start"], all_exons["end"], all_exons["exon_id"]
        ):
            if not transcript_stable_id in result_exons:
                result_exons[transcript_stable_id] = []
                result_exon_ids[transcript_stable_id] = []
            result_exons[transcript_stable_id].append((estart, estop))
            result_exon_ids[transcript_stable_id].append(eid)

        result_exons = pd.Series(
            list(result_exons.values()), index=list(result_exons.keys())
        )
        result_exon_ids = pd.Series(
            list(result_exon_ids.values()), index=list(result_exon_ids.keys())
        )
        result = result.assign(exons=result_exons, exon_stable_ids=result_exon_ids)
        self.sanity_check_transcripts(result)

        return result

    @staticmethod
    def _fix_after_load_df_genes(df):
        if df is not None and len(df) and "transcript_stable_ids" in df.columns:
            res = df.assign(
                transcript_stable_ids=df.transcript_stable_ids.apply(
                    lambda x: tuple(x) if x is not None else None
                ),
            )
            return res
        else:
            return df

    @staticmethod
    def _fix_after_load_df_transcripts(df):
        if df is not None and len(df):
            # assert isinstance(df.exons.iloc[0], np.ndarray)
            # this is what parquet does to our initial tuples. And the downstream expects tuples
            # assert isinstance(df.exons.iloc[0][0], np.ndarray)
            res = df.assign(
                exons=df.exons.apply(lambda x: tuple((tuple(y) for y in x))),
                exon_stable_ids=df.exon_stable_ids.apply(lambda x: tuple(x)),
            )
            assert isinstance(res.exons.iloc[0][0], tuple)
            assert isinstance(res.exons.iloc[0], tuple)
            assert isinstance(res.exon_stable_ids.iloc[0], tuple)
            return res
        else:
            return df

    def sanity_check_transcripts(self, df_transcripts):
        strand_values = set(df_transcripts.strand.unique())
        if strand_values.difference(
            [1, -1]
        ):  # pragma: no cover - defensive, currently handled in gtf parser
            raise ValueError(f"Transcript strand was outside of 1, -1: {strand_values}")

        # can't use self.genes or self.transcript at this point,
        # they rely on df_genes and df_transcripts being set
        genes = df_to_rows(self.df_genes, ["start", "stop"])
        for transcript_stable_id, start, stop, exons, gene_stable_id in zip(
            df_transcripts.index,
            df_transcripts.start,
            df_transcripts.stop,
            df_transcripts.exons,
            df_transcripts.gene_stable_id,
        ):
            if start > stop:
                raise ValueError("start > stop {row}")
            try:
                for estart, estop in exons:
                    if estart < start or estop > stop:
                        raise ValueError(
                            f"Exon outside of transcript: {transcript_stable_id}"
                            f"\ngene was {start}..{stop}"
                            f"\nexon was {estart}..{estop}"
                        )
            except TypeError:  # pragma: no cover
                print(repr((transcript_stable_id, start, stop, exons, gene_stable_id)))
            gene_info = genes[gene_stable_id]
            if start < gene_info.start or stop > gene_info.stop:
                raise ValueError(
                    f"Transcript outside of gene: {transcript_stable_id} {start} {stop} {gene_info.start} {gene_info.stop}"
                )

    def _prepare_df_proteins(self):
        """Get a DataFrame with protein information
        protein_stable_id (index)
        transcript_stable_id,
        gene_stable_id,
        chr,
        strand
        cds - [(start, stop)]  # in genomic coordinates
        """
        gtf = self.get_gtf(["CDS"])
        cds = gtf["CDS"]
        if len(cds) == 0:
            df = pd.DataFrame(
                {
                    "protein_stable_id": [],
                    "transcript_stable_id": [],
                    "gene_stable_id": [],
                    "chr": [],
                    "strand": [],
                    "cds": [],
                }
            ).set_index("protein_stable_id")
            return df
        result = {
            "protein_stable_id": [],
            "transcript_stable_id": [],
            "gene_stable_id": [],
            "chr": [],
            "strand": [],
            "cds": [],
        }

        for protein_stable_id, tuples in dp(cds).groupby("protein_id").itertuples():
            transcript_stable_id = tuples[0].transcript_id
            gene_stable_id = tuples[0].gene_id
            chr = tuples[0].seqname
            strand = tuples[0].strand
            local_cds = list(
                zip((tup.start for tup in tuples), (tup.end for tup in tuples))
            )
            result["protein_stable_id"].append(protein_stable_id[0])
            result["gene_stable_id"].append(gene_stable_id)
            result["transcript_stable_id"].append(transcript_stable_id)
            result["chr"].append(chr)
            result["strand"].append(strand)
            result["cds"].append(local_cds)
        result = pd.DataFrame(result).set_index("protein_stable_id")
        return result

    df_genes = MsgPackProperty(
        lambda self: self.gene_gtf_dependencies,
        lambda self: self.get_additional_gene_gtfs(),
    )
    df_transcripts = MsgPackProperty(
        lambda self: [self.gene_gtf_dependencies, self.job_genes()],
        lambda self: self.get_additional_gene_gtfs(),
    )
    df_proteins = MsgPackProperty(
        lambda self: self.gene_gtf_dependencies,
        lambda self: self.get_additional_gene_gtfs(),
    )

    def get_genes_overlapping(self, chr, start, stop):
        raise ValueError(
            "Use mbf.genomics.Genes.get_overlapping instead. This has no test cases."
        )
        check_overlap = lambda df, interval: np.max(  # noqa: E731
            [
                np.zeros(len(df)),
                np.min(
                    [df.stop.values, np.ones(len(df), dtype=int) * interval[1]], axis=0
                )
                - np.max(
                    [df.start.values, np.ones(len(df), dtype=int) * interval[0]], axis=0
                ),
            ],
            axis=0,
        )
        filter = (self.df_genes["chr"] == chr) & (
            check_overlap(self.df_genes, [start, stop]) > 0
        )
        return self.df_genes[filter]


@class_with_downloads
class HardCodedGenome(GenomeBase):
    def __init__(self, name, chr_lengths, df_genes, df_transcripts, df_proteins):
        super().__init__()
        self.name = name
        self._chr_lengths = chr_lengths
        self._df_genes = GenomeBase._fix_after_load_df_genes(df_genes)
        self._df_transcripts = GenomeBase._fix_after_load_df_transcripts(df_transcripts)
        self._df_proteins = df_proteins

    def get_chromosome_lengths(self):
        return self._chr_lengths.copy()

    def _msg_pack_job(
        self, property_name, filename, callback_function, files_to_invariant_on
    ):
        pass


class GenomePrebuildMixin:
    def _msg_pack_job(
        self, property_name, filename, callback_function, files_to_invariant_on
    ):
        def dump(output_filename):
            df = callback_function(self)
            df.to_parquet(output_filename / filename)

        j = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/{property_name}",
            # we don't use the version for this, since we need it for building
            # various aligner versioned indices
            "4",
            files_to_invariant_on,
            [filename],
            dump,
        )
        j.depends_on_func(property_name, callback_function)
        self._prebuilds.append(j)
        return j

    def build_index(self, aligner, fasta_to_use=None, gtf_to_use=None, postfix=""):
        if fasta_to_use is None:  # pragma: no cover
            _fasta_to_use = "genome.fasta"
            postfix = "genome"
        else:
            _fasta_to_use = fasta_to_use
            if not postfix:
                raise ValueError("Must specify postfix when using fasta_to_use")
        if gtf_to_use is None:  # pragma: no cover
            _gtf_to_use = "genes.gtf"
            gtf_to_use = _gtf_to_use
        else:
            _gtf_to_use = gtf_to_use
            if not postfix:
                raise ValueError("Must specify postfix when using fasta_to_use")
        name = Path(_fasta_to_use).stem

        deps = []
        if hasattr(aligner, "build_index"):
            deps.append(self.find_prebuild(_fasta_to_use))
            if _gtf_to_use:
                deps.append(self.find_prebuild(_gtf_to_use))
            # postfix = ""
            func_deps = {}

            def do_align(output_path):
                aligner.build_index(
                    [self.find_file(_fasta_to_use)],
                    self.find_file(_gtf_to_use) if gtf_to_use is not None else None,
                    output_path,
                )

        elif hasattr(aligner, "build_index_from_genome"):
            if fasta_to_use or gtf_to_use:
                raise ValueError(
                    "Aligner had no build_index, just build_index_from_genome, but fasta_to_use or gtf_to_use were set"
                )
            deps.extend(aligner.get_genome_deps(self))
            func_deps = {
                "build_index_from_genome": aligner.__class__.build_index_from_genome
            }
            postfix += "/" + aligner.get_build_key()

            def do_align(output_path):
                aligner.build_index_from_genome(self, output_path)

        else:
            raise ValueError("Could not find build_index* function")

        min_ver, max_ver = aligner.get_index_version_range()

        aligner_name = aligner.name
        if aligner_name == "STARSolo":
            aligner_name = "STAR"

        job = self.prebuild_manager.prebuild(
            f"{self.prebuild_prefix}/{aligner_name}/{postfix}",
            aligner.version,
            [],
            aligner.get_index_filenames(),
            do_align,
            minimum_acceptable_version=min_ver,
            maximum_acceptable_version=max_ver,
        )
        if hasattr(ppg, "is_ppg2"):
            # that's absolutly a dependency. Good point is though if the output is the same, we will deduplicate
            job.depends_on(
                ppg.ParameterInvariant(job.job_id + "_aligner_ver", aligner.version)
            )
        self.download_genome()  # so that the jobs are there
        job.depends_on(deps)
        for name, f in func_deps.items():
            job.depends_on_func(name, f)
        return job
