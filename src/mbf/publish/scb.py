from pathlib import Path
import bitarray
import struct
import numpy as np
import pypipegraph as ppg
from mbf.externals.util import write_md5_sum
import json
import hashlib


def prep_scb(*objects_to_push, top_level_genes=None):
    """Main entry point to publish data to our scb system"""
    objs = []
    for x in objects_to_push:
        if isinstance(x, dict):
            objs.extend(x.values())
        elif hasattr(x, "__iter__"):
            objs.extend(iter(x))
        elif isinstance(x, File):
            objs.append(x)
        else:
            print(x)
            print(dir(x))
            raise ValueError()
    submission = SCBSubmission(objs, top_level_genes)
    return submission.dump_meta_data()


class _SCBFileEntry:
    def __init__(self, filename_or_job, description, vid, scb_name=None):
        self.description = description
        self.vid = vid
        if not vid:
            raise ValueError("vids are mandatory")
        if (
            isinstance(filename_or_job, tuple)
            and isinstance(filename_or_job[0], ppg.Job)
            and len(filename_or_job) == 2
        ):
            self.path = Path(filename_or_job[1])
            self.deps = filename_or_job[0]
        elif isinstance(filename_or_job, str):
            self.path = Path(filename_or_job)
            self.deps = ppg.FileInvariant(self.path)
        elif isinstance(filename_or_job, Path):
            self.path = filename_or_job
            self.deps = ppg.FileInvariant(self.path)
        elif isinstance(filename_or_job, ppg.Job):
            if hasattr(filename_or_job, "files"):  # ppg2
                fns = filename_or_job.files
            else:
                fns = filename_or_job.filenames  # ppg1
            if len(fns) > 1:
                raise ValueError(
                    "Currently only doing 1 file - pass in a job, filename tuple for teh file you actually want"
                )
            self.path = Path(fns[0])
            self.deps = filename_or_job
        else:
            raise ValueError("Could not handle filename_or_job", repr(filename_or_job))
        self.name = str(self.path)
        if scb_name:
            self.scb_name = scb_name


class File(_SCBFileEntry):
    pass


class BigWigLane(_SCBFileEntry):
    def __init__(self, filename_or_job, genome, vid, scb_name=None):
        super().__init__(filename_or_job, "", vid, scb_name)
        self.genome = genome


class SCBSubmission:
    def __init__(self, objs, top_level_genes=None):
        self.output_path = Path("web/scb")
        self.output_path.mkdir(exist_ok=True, parents=True)
        self.objs = objs
        self.used = set()
        self.genomes = set()
        self.meta_data = {}
        self.deps = []
        self.errors = []
        self.vids = []
        self.master_gene_list_by_genome = {}
        if top_level_genes:
            self.genes = top_level_genes
        else:
            self.genes = {}

        self.genes_to_dump = {}
        self.extract_genomes()
        self.extract_lanes()
        self.extract_genomic_regions()
        self.extract_genes()
        self.extract_single_cell()
        self.extract_tpm_annos()
        self.extract_files()
        if len(self.used) != len(self.objs):
            missing = set(self.objs) - self.used
            raise ValueError("unused objects", missing)
        if self.errors:
            raise ValueError(self.errors)
        self.deps.append(ppg.ParameterInvariant("SCBSubmission_vids", str(self.vids)))

        for x in [
            "extract_genomes",
            "extract_lanes",
            "extract_genomic_regions",
            "extract_genes",
            "extract_single_cell",
            "extract_tpm_annos",
            "extract_files",
        ]:
            self.deps.append(
                ppg.FunctionInvariant(f"SCBSubmission.{x}", getattr(self.__class__, x))
            )

    def dump_meta_data(self):
        self.dump_meta_data_json()
        self.dump_rsync_list()

    def dump_meta_data_json(self):
        output_filename = "web/scb/metadata.json"

        def dump(output_filename):
            descend_and_replace_callbacks(self.meta_data)
            with open(output_filename, "w") as op:
                json.dump(
                    self.meta_data, op, indent=4
                )  # make sure it's sane, no funny objects

        return ppg.FileGeneratingJob(output_filename, dump).depends_on(self.deps)

    def dump_rsync_list(self):
        output_filename = "web/scb/rsync_list.txt"

        def dump(output_filename):
            paths_to_copy = set()
            for group in self.meta_data:
                for entry in self.meta_data[group]:
                    for x in [  # this is bad style.
                        "table_path",
                        "path_bigbed",
                        "path_table",
                        "path_bam",
                        "path_bw",
                        "path",
                        "path_json",
                        "path_json_idx",
                    ]:
                        if x in entry:  # genes
                            paths_to_copy.add(Path(entry[x]).absolute())
                    if "path_bam" in entry:
                        paths_to_copy.add(Path(entry["path_bam"] + ".bai").absolute())
                    if 'paths' in entry:
                        for p in entry['paths']:
                            paths_to_copy.add(Path(p).absolute())

            for fn in Path("web/scb").glob("*"):
                paths_to_copy.add(fn.absolute())
            output = ""
            # paths_to_copy.add("/project/web/scb/metadata.json")
            for x in sorted([str(x) for x in paths_to_copy]):
                if x.startswith("/project"):
                    output += x[len("/project/") :]
                    output += "\n"
            Path(output_filename).write_text(output)

        return (
            ppg.FileGeneratingJob(output_filename, dump)
            .depends_on(self.deps)
            .depends_on(self.dump_meta_data_json())
        )

    def extract_genomes(self):
        for o in self.objs:
            if hasattr(o, "genome"):
                self.genomes.add(o)

    def extract_lanes(self):
        self.meta_data["lanes"] = []
        self.meta_data["lanes_bigwig"] = []
        for entry in self.objs:
            if hasattr(entry, "get_bam_names"):
                vids = flatten_vid(entry.vid, entry, self.errors)
                self.vids.append(vids)
                fn = entry.get_bam_names()[0]
                md5 = get_md5(fn, entry.load())
                self.meta_data["lanes"].append(
                    {
                        "vids": vids,
                        "path_bam": fn,
                        "displayname": get_scb_name(entry),
                        "md5": lambda md5=md5: Path(md5[1]).read_text(),
                        "md5_path": fn,
                        "scb_comment": (
                            entry.scb_comment if hasattr(entry, "scb_comment") else ""
                        ),
                        "id": (
                            entry.scb_id if hasattr(entry, "scb_id") else None
                        ),  # this allows you to permanently tie it to an entry in SCBs databases
                        # "browser_options": entry.gbrowse_options,
                        # "has_splicing": entry.has_spliced_reads,
                        "genome_label": extract_genome_label(entry),
                    }
                )
                self.deps.append(md5[0])
                self.register_used(entry)
            elif isinstance(entry, BigWigLane):
                vids = flatten_vid(entry.vid, entry, self.errors)
                self.vids.append(vids)
                fn = str(entry.path.absolute().relative_to("/project"))
                md5 = get_md5(fn, entry.deps)
                self.meta_data["lanes_bigwig"].append(
                    {
                        "vids": vids,
                        "path_bw": fn,
                        "displayname": get_scb_name(entry),
                        "md5": lambda md5=md5: Path(md5[1]).read_text(),
                        "md5_path": fn,
                        "scb_comment": (
                            entry.scb_comment if hasattr(entry, "scb_comment") else ""
                        ),
                        "id": (
                            entry.scb_id if hasattr(entry, "scb_id") else None
                        ),  # this allows you to permanently tie it to an entry in SCBs databases
                        "genome_label": extract_genome_label(entry),
                    }
                )
                self.deps.append(md5[0])
                self.register_used(entry)

    def extract_tpm_annos(self):
        self.meta_data["tag_count_annos"] = []
        for entry in self.objs:
            for kind, class_name in [
                ("TPM", "NormalizationTPM"),
                ("CPM", "NormalizationCPM"),
            ]:
                if (
                    entry.__class__.__name__ == class_name
                ):  # don't want to import this here
                    vids = flatten_vid(entry.vid, entry, self.errors)
                    self.vids.append(vids)

                    def calc_md5(entry=entry):
                        genes = self.get_genes(entry.genome)
                        return hashlib.md5(
                            (
                                "".join([str(x) for x in genes.df[entry.columns[0]]])
                            ).encode("utf-8")
                        ).hexdigest()

                    job_table = self.table_jobs[entry.genome]
                    fn_table = Path(job_table.filenames[0])

                    self.meta_data["tag_count_annos"].append(
                        {
                            "vids": vids,
                            "displayname": get_scb_name(
                                entry
                            ),  # currently, these are not being shown in the SCB, we use the sample info instead...
                            "column_name": entry.columns[0],
                            "md5": calc_md5,
                            "table_path": str(
                                fn_table.absolute().relative_to("/project")
                            ),
                            "scb_comment": getattr(entry, "scb_comment", ""),
                            "kind": kind,
                            "id": getattr(
                                entry, "scb_id", None
                            ),  # this allows you to permanently tie it to an entry in SCBs databases
                            "genome_label": extract_genome_label(entry),
                        }
                    )
                    genes = self.get_genes(entry.genome)
                    if not genes.has_annotator(entry):
                        raise ValueError(
                            f"Genes object {genes.name} did not have annotator {entry.columns[0]}."
                            "\nprep_scb won't add it!"
                            "\nPerhaps you need to pass in prep_scb(..., top_level_genes={genome: genes}) to use?"
                            "\n(e.g. if you have added tpm annos just to protein_coding genes)"
                        )
                    self.deps.append(genes.add_annotator(entry))
                    self.register_used(entry)

    def extract_genomic_regions(self):
        self.meta_data["genomic_regions"] = []
        for entry in self.objs:
            if entry.__class__.__name__ == "GenomicRegions":
                vids = flatten_vid(entry.vid, entry, self.errors)
                self.vids.append(vids)
                job_bigbed, fn_bigbed = entry.write_bigbed()
                job_table, fn_table = entry.write()
                md5 = get_md5(fn_bigbed, job_bigbed)
                self.meta_data["genomic_regions"].append(
                    {
                        "vids": vids,
                        "path_bigbed": str(
                            fn_bigbed.absolute().relative_to("/project")
                        ),
                        "path_table": str(fn_table.absolute().relative_to("/project")),
                        "displayname": get_scb_name(entry),
                        "md5": lambda md5=md5: Path(md5[1]).read_text(),
                        "md5_path": str(
                            fn_bigbed.absolute().relative_to("/project")
                        ),  # this is a bit more robust than taking the table path, since it does not change just because columns were added
                        "scb_comment": entry.scb_comment
                        if hasattr(entry, "scb_comment")
                        else "",
                        "id": entry.scb_id
                        if hasattr(entry, "scb_id")
                        else None,  # this allows you to permanently tie it to an entry in SCBs databases
                        "size": lambda entry=entry: len(entry.df),
                        "genome_label": extract_genome_label(entry),
                    }
                )
                self.deps.append(entry.load)
                self.deps.append(job_bigbed)
                self.deps.append(job_table)
                self.deps.append(md5[0])
                self.register_used(entry)

    def extract_genes(self):
        self.meta_data["genes"] = []
        self.table_jobs = {}
        # collect first, so that write_gene_sql can make
        # a parameter invariant over all of them.
        genomes = set()
        for entry in self.objs:
            if entry.__class__.__name__ == "Genes":
                if not entry.genome in self.genes_to_dump:
                    self.genes_to_dump[entry.genome] = []
                    genomes.add(entry.genome)
                self.genes_to_dump[entry.genome].append(entry)

        for g in genomes:
            self.table_jobs[g] = self.write_gene_sql(g)

        for entry in self.objs:
            if entry.__class__.__name__ == "Genes":
                vids = flatten_vid(entry.vid, entry, self.errors)
                self.vids.append(vids)

                def calc_md5(entry=entry):
                    # parent row is the row(=index) in the very top level parent!
                    # but we need a bool vector wheter it's
                    h = hashlib.md5()
                    for stable_id in sorted(entry.df.gene_stable_id):
                        h.update(stable_id.encode("utf-8"))
                    return h.hexdigest()

                scb_comment = getattr(entry, "scb_comment", "")
                if scb_comment is None:
                    scb_comment = ""
                fn_table = Path(self.table_jobs[entry.genome].filenames[0])

                self.meta_data["genes"].append(
                    {
                        "vids": vids,
                        "displayname": get_scb_name(entry),
                        "dbname": entry.name,
                        "md5": calc_md5,
                        "scb_comment": scb_comment,
                        "id": getattr(
                            entry, "scb_id", None
                        ),  # this allows you to permanently tie it to an entry in SCBs databases
                        "table_path": str(fn_table.absolute().relative_to("/project")),
                        "genome_label": extract_genome_label(entry),
                    }
                )
                self.deps.append(entry.load())
                self.deps.append(self.get_genes(entry.genome).load())
                self.deps.append(self.table_jobs[entry.genome])
                self.register_used(entry)

    def extract_single_cell(self):
        self.meta_data["single_cell_cirrocumulus"] = []
        for entry in self.objs:
            if entry.__class__.__name__ == "SingleCellForSCB":
                vids = flatten_vid(entry.vid, entry, self.errors)
                self.vids.append(vids)

                md5 = get_md5(entry.output_filenames[0], entry.load())
                self.meta_data["single_cell_cirrocumulus"].append(
                    {
                        "vids": vids,
                        "path_h5ad": str(entry.h5ad),
                        "path_json": str(entry.output_filenames[0]),
                        "path_json_idx": str(entry.output_filenames[1]),
                        "paths": [str(x) for x in entry.output_filenames[2:]],
                        "displayname": get_scb_name(entry),
                        "md5": lambda md5=md5: Path(md5[1]).read_text(),
                        "md5_path": str(
                            entry.output_filenames[0].absolute().relative_to("/project")
                        ),
                        "scb_comment": entry.scb_comment
                        if hasattr(entry, "scb_comment")
                        else "",
                        "id": entry.scb_id
                        if hasattr(entry, "scb_id")
                        else None,  # this allows you to permanently tie it to an entry in SCBs databases
                        "genome_label": extract_genome_label(entry),
                        "stats": entry.stats,
                    }
                )
                self.deps.append(entry.load)
                self.deps.append(md5[0])
                self.register_used(entry)

    def extract_files(self):
        self.meta_data["files"] = []
        for entry in self.objs:
            if entry.__class__ == File:
                vids = flatten_vid(entry.vid, entry, self.errors)
                self.vids.append(vids)

                md5 = get_md5(entry.path, entry.deps)  # todo : grab from ppg2?
                rel_path = str(entry.path.absolute().relative_to("/project"))
                self.meta_data["files"].append(
                    {
                        "vids": vids,
                        "path": rel_path,
                        "displayname": get_scb_name(entry),
                        "md5": lambda md5=md5: Path(md5[1]).read_text(),
                        "md5_path": rel_path,
                        "description": entry.description,
                        "scb_comment": entry.scb_comment
                        if hasattr(entry, "scb_comment")
                        else "",
                        "id": entry.scb_id
                        if hasattr(entry, "scb_id")
                        else None,  # this allows you to permanently tie it to an entry in SCBs databases
                    }
                )
                self.deps.append(entry.deps)
                self.deps.append(md5[0])  # that's a job :)
                self.register_used(entry)

    def register_used(self, entry):
        if entry in self.used:
            raise ValueError("double use", entry)
        self.used.add(entry)

    def get_genes(self, genome):

        if not genome in self.genes:
            import mbf.genomics

            self.genes[genome] = mbf.genomics.genes.Genes(genome)
        return self.genes[genome]

    def write_gene_sql(self, genome):  # noqa:C901
        import sqlite3

        genes = self.get_genes(genome)

        output_filename = Path("web", "scb", "genes", genes.name + ".sqlite_split_df")
        output_filename.parent.mkdir(exist_ok=True, parents=True)

        def chunks(ll, n):
            n = max(1, n)
            return [ll[i : i + n] for i in range(0, len(ll), n)]

        def write(output_filename, genes=genes):
            # column_properties = self.get_column_properties()
            df = genes.df  # in the order defined by priority

            for col in df.columns:
                if df[col].dtype == "object" or col in ["chr", "biotype"]:
                    # print("fixing", col)
                    df = df.assign(**{col: [str(x) for x in df[col].values]})
            df = df.rename(columns={"gene_stable_id": "stable_id"})
            # print(df.columns)

            index_columns = ["stable_id"]
            output = {}
            index_arg = {}
            column_to_df_no = {c: 0 for c in index_columns}
            for ii, columns in enumerate(chunks(list(df.columns), 990)):
                for idx_column in index_columns:
                    if idx_column not in columns:
                        columns.append(idx_column)
                sub_df = df[columns]
                output["df_%i" % ii] = sub_df
                index_arg["df_%i" % ii] = index_columns
                for c in columns:
                    if c not in index_columns:
                        column_to_df_no[c] = ii
            # if Path(output_filename).exists():
            # print("unlinking", output_filename)
            # Path(output_filename).unlink()
            with sqlite3.connect(output_filename, isolation_level="DEFERRED") as conn:
                cur = conn.cursor()
                cur.execute("BEGIN TRANSACTION")
                for table_name, sub_df in output.items():
                    # print(sub_df.dtypes)
                    sub_df.to_sql(table_name, conn)
                    statement = (
                        "CREATE TABLE '%s' (column_name TEXT, column_type TEXT)"
                        % (table_name + "_meta",)
                    )
                    cur.execute(statement)
                    meta = []
                    for column_name, column_type in sub_df.dtypes.iteritems():
                        meta.append((column_name, str(column_type)))
                    cur.executemany(
                        "INSERT into '%s_meta' VALUES (?, ?)" % (table_name,), meta
                    )
                    indexed_columns = index_arg[table_name]
                    if indexed_columns:  # index columns are in every table
                        for column_name in indexed_columns:
                            # if column_properties[column_name].get("nocase", False):
                            index_options = "collate nocase"
                            # else:
                            #    index_options = ""
                            cur.execute(
                                "CREATE INDEX '%s_%s' ON '%s' ('%s' %s)"
                                % (
                                    table_name,
                                    column_name,
                                    table_name,
                                    column_name,
                                    index_options,
                                )
                            )
                conn.commit()
                cur.execute("BEGIN TRANSACTION")

                # store the mapping column_name -> df no
                cur.execute(
                    """CREATE TABLE column_to_df (
                            column_name VARCHAR,
                            df_no integer,
                            overall_column_pos integer
                )"""
                )
                d = []
                col_list = list(df.columns)
                for column_name, df_no in column_to_df_no.items():
                    d.append((column_name, df_no, col_list.index(column_name)))
                cur.executemany(
                    "INSERT into column_to_df (column_name, df_no, overall_column_pos) values (?, ?, ?)",
                    d,
                )

                # now store the column_properties

                cur.execute(
                    """CREATE TABLE column_properties (
                            column_name VARCHAR,
                            property_name VARCHAR,
                            value VARCHAR
                )"""
                )
                column_properties = {
                    "chr": {
                        "user_visible": False,
                        "priority": -997,
                        "description": "On which chromosome (or contig) the gene is loacted",
                    },
                    "start": {
                        "user_visible": False,
                        "priority": -997,
                        "description": "Left most position of this gene",
                    },
                    "stop": {
                        "user_visible": False,
                        "priority": -997,
                        "description": "Right most position of this gene",
                    },
                    "tss": {
                        "user_visible": False,
                        "priority": -997,
                        "description": "Position of Transcription start site on the chromosome",
                    },
                    "tes": {
                        "user_visible": False,
                        "priority": -997,
                        "description": "Position of the annotated end of transcription",
                    },
                    "stable_id": {
                        "user_visible": True,
                        "priority": -1000,
                        "index": True,
                        "description": "Unique identification name for this gene",
                    },
                    "name": {
                        "user_visible": True,
                        "priority": -999,
                        "index": True,
                        "description": "Offical name for this gene",
                        "nocase": True,
                    },
                    "strand": {
                        "user_visible": True,
                        "priority": -998,
                        "description": "Which is the coding strand (1 forward, -1 reverse)",
                    },
                    "description": {
                        "user_visible": True,
                        "priority": -998.5,
                        "description": "A short description of the gene",
                    },
                    "transcript_stable_ids": {
                        "user_visible": False,
                        "priority": 999,
                        "description": "Transcripts of this gene",
                    },
                }
                column_properties = {
                    k: v for (k, v) in column_properties.items() if k in df.columns
                }
                for k in df.columns:
                    if not k in column_properties:
                        column_properties[k] = {"priority": 0}
                for (col, properties) in column_properties.items():
                    for prop_name, prop_value in properties.items():
                        cur.execute(
                            "INSERT into column_properties VALUES (?, ?, ?)",
                            (col, prop_name, prop_value),
                        )
                conn.commit()

                valid_types = (
                    1,  # sign, up=1, down=-1, both=0
                    2,  # abs(value) >= x
                    3,  # filter to text value
                    4,  # minimum(columns) >= x
                    5,  # x <= value
                )
                conn.text_factory = str
                cur = conn.cursor()
                cur.execute("BEGIN TRANSACTION")
                cur.execute("DROP TABLE IF EXISTS subsets")
                cur.execute(
                    """CREATE TABLE subsets (subset_name text, sheet_name text,  included_vector blob)"""
                )
                cur.execute("DROP TABLE IF EXISTS subset_filter_columns")
                cur.execute(
                    """CREATE TABLE subset_filter_columns (subset_name text, column_desc text, column_name text, filter_type int, filter_default_value text )"""
                )
                cur.execute("DROP TABLE IF EXISTS subset_relevant_columns")
                cur.execute(
                    """CREATE TABLE subset_relevant_columns (subset_name text, column_name text)"""
                )
                stable_id_to_no = {}
                for ii, stable_id in enumerate(
                    sorted(genes.df["gene_stable_id"])
                ):  # the vectors are defined as to be on the sorted list of genes
                    stable_id_to_no[stable_id] = ii

                for entry in self.genes_to_dump[genome]:
                    if entry.genome != genome:
                        raise ValueError("or this is it")
                    if entry.genome != genes.genome:
                        raise ValueError("or this is it2")
                    entry_vector = np.zeros((len(genes.df)), dtype=bool)
                    for stable_id in entry.df["gene_stable_id"]:
                        entry_vector[stable_id_to_no[stable_id]] = 1

                    ba = bitarray.bitarray()
                    ba.extend(entry_vector)
                    vector_string = (
                        b"B" + struct.pack(b"I", len(entry_vector)) + ba.tobytes()
                    )  # vector_string = unicode("".join(['1' if x else '0' for x in entry_vector]))
                    # if len(vector_string) - 1!= len(self.master_genes.df):
                    # print len(vector_string), len(self.master_genes.df)
                    # raise ValueError("Entry vector length unequal gene list length")

                    cur.execute(
                        "INSERT into subsets (subset_name, sheet_name, included_vector) VALUES (?, ?, ?)",
                        (
                            entry.name,
                            entry.sheet_name if entry.sheet_name else "default",
                            memoryview(vector_string),
                        ),
                    )
                    if hasattr(entry, "further_filter_columns"):
                        for (
                            column_desc,
                            column_name,
                            column_type,
                            default_value,
                        ) in entry.further_filter_columns:
                            if not column_type in valid_types:
                                raise ValueError(
                                    "Invalid column type: %s" % column_type
                                )
                            if column_type == 4:  # stupid fileformats are stupid
                                column_name = json.dumps(column_name)
                            print(
                                "dbg",
                                entry.name,
                                column_desc,
                                column_name,
                                column_type,
                                default_value,
                            )
                            print(
                                "type",
                                type(entry.name),
                                type(column_desc),
                                type(column_name),
                                type(column_type),
                                type(default_value),
                            )

                            cur.execute(
                                "INSERT into subset_filter_columns VALUES (?, ?, ?, ?, ?)",
                                (
                                    entry.name,
                                    column_desc,
                                    column_name,
                                    column_type,
                                    default_value,
                                ),
                            )
                    if hasattr(entry, "subset_relevant_columns"):
                        cur.executemany(
                            "INSERT into subset_relevant_columns VALUES (?,?)",
                            (
                                [
                                    (entry.name, column_name)
                                    for column_name in entry.subset_relevant_columns
                                ]
                            ),
                        )
                conn.commit()

        subset_and_filter_columns = {}
        subset_and_filter_columns[genome.name] = {}
        for x in self.genes_to_dump[genome]:
            subset_and_filter_columns[x.name] = (
                getattr(x, "subset_relevant_columns", []),
                getattr(x, "further_filter_columns", []),
            )

        return (
            ppg.FileGeneratingJob(output_filename, write)
            .depends_on(
                genes.annotate(), [x.load() for x in self.genes_to_dump[genome]]
            )
            .depends_on(
                ppg.ParameterInvariant(
                    output_filename, ppg.util.freeze(subset_and_filter_columns)
                )
            )
        )


def get_md5(input_name, input_job):
    input_name = Path(input_name)
    target_name = input_name.with_name(input_name.name + ".md5sum")

    def gen(target_filename):
        target_filename = Path(target_filename)
        write_md5_sum(target_filename.with_name(target_filename.stem))

    return (
        ppg.FileGeneratingJob(target_name, gen).depends_on(input_job),
        str(target_name),
    )


def flatten_vid(vid, entry, errors, collector=None):
    if collector is None:
        collector = []
    if vid is None or vid == []:
        errors.append(
            (
                entry.name if hasattr(entry, "name") else entry.column_name,
                "Vid was none - all samples have to be associated with a vid",
            )
        )
    elif hasattr(vid, "__iter__") and not isinstance(vid, str):
        for v in vid:
            try:
                flatten_vid(v, entry, errors, collector)
            except RecursionError:
                print(entry)
                raise
    else:
        collector.append(vid)
    return collector


def get_scb_name(obj):
    if hasattr(obj, "scb_name") and obj.scb_name:
        return obj.scb_name
    elif hasattr(obj, "sheet_name") and obj.sheet_name:
        return obj.sheet_name + "/" + obj.name
    elif hasattr(obj, "plot_name") and obj.plot_name:
        return obj.plot_name
    else:
        return obj.name


def extract_genome_label(obj):
    if obj.genome.species == "Mus_musculus":
        return "mm10"
    elif obj.genome.species == "Homo_sapiens":
        if int(obj.genome.revision) >= 77:
            return "hg38"
        else:
            return "hg19"
    elif obj.genome.species == "Drosophila_melanogaster":
        return "dm89"
    elif obj.genome.species == "Ustilago_maydis":
        return "um33"
    else:
        raise ValueError(
            "don't know how to tell the scb about this genome: %s" % obj.genome
        )


def descend_and_replace_callbacks(d):
    for k, v in d.items():
        if hasattr(v, "__call__"):
            d[k] = v()
        elif isinstance(v, dict):
            descend_and_replace_callbacks(v)
        elif isinstance(v, list):
            for vx in v:
                if isinstance(vx, dict):
                    descend_and_replace_callbacks(vx)
