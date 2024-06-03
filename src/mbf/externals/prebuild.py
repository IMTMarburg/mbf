"""
Many algorithms need prebuild data structures (indices and so on)
which both are too time consuming to build, to big to copy to each
experiment and need to be versioned,
but they can often be shared among versions.

This has drastically simplified with ppg2.


"""

import socket
from .util import Version, sort_versions, UpstreamChangedError, write_md5_sum
import pypipegraph as ppg
from pathlib import Path
import time
import stat
import os
import json
import contextlib


if hasattr(ppg, "is_ppg2"):  # noqa: C901
    import pypipegraph2 as ppg2

    PrebuildFunctionInvariantFileStoredExploding = ppg.FunctionInvariant
    _PrebuildFileInvariantsExploding = lambda job_id, filenames: [  # noqa:E731
        ppg.FileInvariant(x) for x in filenames
    ]
    PrebuildJob = ppg2.SharedMultiFileGeneratingJob
else:

    class PrebuildFunctionInvariantFileStoredExploding(ppg.FunctionInvariant):
        def __init__(self, storage_filename, func):
            self.is_prebuild = True
            super().__init__(storage_filename, func)

        @classmethod
        def hash_function(cls, function):
            new_source, new_funchash, new_closure = cls._hash_function(function)
            return cls._compare_new_and_old(
                new_source, new_funchash, new_closure, False
            )

        def _get_invariant(self, old, all_invariant_stati):
            stf = Path(
                self.job_id
            )  # the old file format - using just the function's dis-ed code.
            stf2 = Path(self.job_id).with_name(
                stf.name + "2"
            )  # the new style, dict based storage just like FunctionInvariant after 0.190
            new_source, new_func_hash, new_closure = self._hash_function(self.function)
            if stf2.exists():
                old_hash = json.loads(stf2.read_text())
            elif stf.exists():
                old_hash = stf.read_text()
                new_closure = ""
            else:
                new_value = self._compare_new_and_old(
                    new_source, new_func_hash, new_closure, False
                )
                stf2.write_text(json.dumps(new_value))
                return old  # signal no change necessary.

            try:
                new_hash = self._compare_new_and_old(
                    new_source, new_func_hash, new_closure, old_hash
                )
                if new_hash != old_hash:
                    self.complain_about_hash_changes(new_hash)
                else:
                    return old
            except ppg.NothingChanged as e:
                # we accept the stuff there as no change.
                # and we write out the new value, because it might be a format change.
                try:
                    stf2.write_text(json.dumps(e.new_value))
                except OSError as e2:
                    if "Read-only file system" in str(e2):
                        import warnings

                        warnings.warn(
                            "PrebuildFunctionInvariantFileStoredExploding: Could not update %s to newest version - read only file system"
                            % stf
                        )
                raise e
            raise NotImplementedError("Should not happen")

        def complain_about_hash_changes(self, invariant_hash):
            stf = Path(self.job_id)
            try:
                of = stf.with_name(stf.name + ".changed")
                of.write_text(json.dumps(invariant_hash))
            except IOError:  # noqa: E722 pragma: no cover
                # fallback if the stf directory is not writeable.
                of = Path(stf.name + ".changed")  # pragma: no cover
                of.write_text(json.dumps(invariant_hash))  # pragma: no cover
            raise UpstreamChangedError(
                (
                    "Calculating function changed.\n"
                    "If you are actively working on it, you need to bump the version:\n"
                    "If not, you need to figure out what's causing the change.\n"
                    "Do not nuke the job info (%s) light heartedly\n"
                    "To compare, run \n"
                    "icdiff %s %s"
                )
                % (self.job_id, Path(self.job_id).absolute(), of.absolute())
            )

    class _PrebuildFileInvariantsExploding(ppg.MultiFileInvariant):
        """Used by PrebuildJob to handle input file deps"""

        def __new__(cls, job_id, filenames):
            job_id = "PFIE_" + str(job_id)
            return ppg.Job.__new__(cls, job_id)

        def __init__(self, job_id, filenames):
            job_id = "PFIE_" + str(job_id)
            self.filenames = filenames
            for f in filenames:
                if not (isinstance(f, str) or isinstance(f, Path)):  # pragma: no cover
                    raise ValueError(f"filenames must be str/path. Was {repr(f)}")
            self.is_prebuild = True
            ppg.Job.__init__(self, job_id)

        def calc_checksums(self, old):
            """return a list of tuples
            (filename, filetime, filesize, checksum)"""
            result = []
            if old:
                old_d = {x[0]: x[1:] for x in old}
            else:
                old_d = {}
            for fn in self.filenames:
                if not os.path.exists(fn):
                    result.append((fn, None, None, None))
                else:
                    st = os.stat(fn)
                    filetime = st[stat.ST_MTIME]
                    filesize = st[stat.ST_SIZE]
                    if (
                        fn in old_d
                        and (old_d[fn][0] == filetime)
                        and (old_d[fn][1] == filesize)
                    ):  # we can reuse the checksum
                        result.append((fn, filetime, filesize, old_d[fn][2]))
                    else:
                        result.append(
                            (fn, filetime, filesize, ppg.util.checksum_file(fn))
                        )
            return result

        def _get_invariant(self, old, all_invariant_stati):
            if not old:
                old = self.find_matching_renamed(all_invariant_stati)
            checksums = self.calc_checksums(old)
            if old is False:
                raise ppg.ppg_exceptions.NothingChanged(checksums)
            # elif old is None: # not sure when this would ever happen
            # return checksums
            else:
                old_d = {x[0]: x[1:] for x in old}
                checksums_d = {x[0]: x[1:] for x in checksums}
                for fn in self.filenames:
                    if (
                        fn in old_d
                        and old_d[fn][2] != checksums_d[fn][2]
                        and old_d[fn][2] is not None
                    ):
                        raise UpstreamChangedError(
                            """Upstream file changed for job, bump version or rollback.
Job: %s
File: %s"""
                            % (self, fn)
                        )
                    # return checksums
            raise ppg.ppg_exceptions.NothingChanged(checksums)

    class PrebuildJob(ppg.MultiFileGeneratingJob):
        def __new__(cls, filenames, calc_function, output_path):
            if not hasattr(filenames, "__iter__"):
                raise TypeError("filenames was not iterable")
            for x in filenames:
                if not (isinstance(x, str) or isinstance(x, Path)):
                    raise TypeError(
                        "filenames must be a list of strings or pathlib.Path"
                    )
            for of in filenames:
                if Path(of).is_absolute():
                    raise ValueError("output_files must be relative")
            filenames = cls._normalize_output_files(filenames, output_path)
            job_id = ":".join(sorted(str(x) for x in filenames))
            res = ppg.Job.__new__(cls, job_id)
            res.filenames = filenames
            res.output_path = Path(output_path)
            return res

        @classmethod
        def _normalize_output_files(cls, output_files, output_path):
            output_files = [
                Path(cls.verify_job_id(output_path / of)) for of in output_files
            ]
            output_files.append(Path(cls.verify_job_id(output_path / "mbf.done")))
            return output_files

        def __init__(self, output_files, calc_function, output_path):
            output_files = self._normalize_output_files(output_files, output_path)
            output_path.mkdir(parents=True, exist_ok=True)

            self.real_callback = calc_function
            self.is_prebuild = True

            def calc():
                self.real_callback(output_path)
                output_files[-1].write_text(str(time.time()))
                for fn in output_files[:-1]:
                    if os.path.exists(fn):
                        write_md5_sum(fn)

            super().__init__(output_files, calc, rename_broken=True, empty_ok=True)
            self.output_path = output_path

        def depends_on_func(self, name, func):
            job = PrebuildFunctionInvariantFileStoredExploding(
                self.output_path / ("%s.md5sum" % (name,)), func
            )
            self.depends_on(job)
            return job

        def depends_on_file(self, filename):
            job = _PrebuildFileInvariantsExploding(filename, [filename])
            self.depends_on(job)
            return job

        def depends_on(self, jobs):
            for job in ppg.util.flatten_jobs(jobs):
                if not hasattr(job, "is_prebuild") or not job.is_prebuild:
                    raise ppg.JobContractError(
                        "%s depended on a non-prebuild dependency %s - not supported"
                        % (self, job)
                    )
                ppg.Job.depends_on(self, job)
            return self

        def inject_auto_invariants(self):
            self.depends_on_func("mbf.func", self.real_callback)

        def invalidated(self, reason):
            exists = [Path(of).exists() for of in self.filenames]
            if all(exists) or not any(exists):
                pass
            else:
                raise ValueError(
                    "Some output files existed, some don't - undefined state, manual cleanup needed\n:%s"
                    % (list(zip(self.filenames, exists)))
                )
            self.was_invalidated = True

        def name_file(self, output_filename):
            """Adjust path of output_filename by job path"""
            return self.output_path / output_filename

        def find_file(self, output_filename):
            """Search for a file named output_filename in the job's known created files"""
            of = self.name_file(output_filename)
            for fn in self.filenames:
                if of.resolve() == Path(fn).resolve():
                    return of
            else:
                raise KeyError("file not found: %s" % output_filename)

        def __getitem__(self, item):
            return self.find_file(item)


class DummyJob:
    """just enough of the Jobs interface to ignore the various calls
    and allow finding the msgpack jobs
    """

    def __init__(self, output_path, filenames):
        self.output_path = output_path
        self.filenames = filenames
        # self.job_id = ":".join(sorted(str(x) for x in filenames))

    def __str__(self):
        return "mbf.externals.prebuild.DummyJob(%s) %s" % (
            self.output_path,
            self.filenames,
        )

    def depends_on(self, _other_job):  # pragma: no cover
        return self

    def depends_on_func(self, _name, _func):  # pragma: no cover
        return self

    def depends_on_file(self, _filename):  # pragma: no cover
        return self

    def depends_on_params(self, _values):  # pragma: no cover
        return self

    def name_file(self, output_filename):
        """Adjust path of output_filename by job path"""
        return self.output_path / output_filename

    def find_file(self, output_filename):
        """Search for a file named output_filename in the job's known created files"""
        of = self.name_file(output_filename)
        for fn in self.filenames:
            if of.resolve() == Path(fn).resolve():
                return of
        else:
            raise KeyError(
                f"file not found: {output_filename}, searching in {self.output_path}"
            )

    def __iter__(self):
        yield self


class PrebuildManager:
    def __init__(self, prebuilt_path, hostname=None):
        self._prebuilt_path = Path(prebuilt_path)
        self.hostname = hostname if hostname else socket.gethostname()
        (self.prebuilt_path / self.hostname).mkdir(exist_ok=True, parents=True)

    @property
    def prebuilt_path(self):
        if hasattr(ppg, "is_ppg2"):
            return self._prebuilt_path / "ppg2"
        else:
            return self._prebuilt_path

    def _find_versions(self, name):
        result = {}
        dirs_to_consider = [
            p
            for p in self.prebuilt_path.glob("*")
            if (p / name).exists() and p.name != self.hostname
        ]
        # prefer versions from this host - must be last!
        dirs_to_consider.append(self.prebuilt_path / self.hostname)
        for p in dirs_to_consider:
            for v in (p / name).glob("*"):
                if (v / "mbf.done").exists():
                    result[v.name] = v
        return result

    def prebuild(  # noqa: C901
        self,
        name,
        version,
        input_files,
        output_files,
        calculating_function,
        minimum_acceptable_version=None,
        maximum_acceptable_version=None,
        further_function_deps={},
        ppg2_resources=None,
        remove_unused=True,
    ):
        """Create a job that will prebuilt the files if necessary

        @further_function_deps is a dictionary name => func,

        PPG1
            : version, minimum_acceptable_version, maximum_acceptable_version
            will be used to find a compatible version (if possible).

            further_function_deps will end up as PrebuildFunctionInvariantFileStoredExploding
            in the correct directory.

            ppg2_resources and remove_unused is ignored

        PPG2:
            version is a ParameterInvariant.
            If the version was never build ppg2.SharedMultiFileGeneratingJob
            will be used to build it. That is then supposed to map it to
            a symlink, if the output is identical.

            That means every used version will be built at least once.
            But you don't have to think about compability.

            minimum_acceptable_version and maximum_acceptable_version are ignored in that contexxt

            remove_unused is passed to SharedMultiFileGeneratingJob


        If used outside of ppg, it will not build anything, complain if it can't find it,
        and give something back that's just about ok to use in mbf.genomes msgpack jobs.

        todo: cores needed for ppg2

        """
        if hasattr(ppg, "is_ppg2"):
            import pypipegraph2 as ppg2

            if ppg2_resources is None:
                ppg2_resources = ppg2.Resources.SingleCore
            return self._prebuild_ppg2(
                name,
                version,
                input_files,
                output_files,
                calculating_function,
                further_function_deps,
                ppg2_resources,
                remove_unused,
            )
        else:
            return self._prebuild_ppg1(
                name,
                version,
                input_files,
                output_files,
                calculating_function,
                minimum_acceptable_version,
                maximum_acceptable_version,
                further_function_deps,
            )

    def _prebuild_ppg2(
        self,
        name,
        version,
        input_files,
        output_files,
        calculating_function,
        further_function_deps,
        resources,
        remove_unused,
    ):
        # import pypipegraph2 as ppg2

        if isinstance(output_files, (str, Path)):
            output_files = [output_files]
        output_files = [Path(of) for of in output_files]
        output_dir_prefix = self.prebuilt_path / self.hostname / name
        if ppg.inside_ppg():
            return self._prebuild_ppg2_ppg(
                name,
                version,
                input_files,
                output_files,
                calculating_function,
                further_function_deps,
                resources,
                remove_unused,
                output_dir_prefix,
            )
        else:
            return self._prebuild_ppg2_no_ppg(
                name,
                version,
                input_files,
                output_files,
                calculating_function,
                further_function_deps,
                resources,
                remove_unused,
                output_dir_prefix,
            )

    def _prebuild_ppg2_ppg(
        self,
        name,
        version,
        input_files,
        output_files,
        calculating_function,
        further_function_deps,
        resources,
        remove_unused,
        output_dir_prefix,
    ):
        import pypipegraph2 as ppg2

        def adapt_calling(output_files, output_dir_prefix):
            # the old just took a prefix to put it's files in,
            # and the file names had to be captured from the outside.
            # the new one, a bit more sensibly, I suppose, receives both.
            return calculating_function(output_dir_prefix)

        adapt_calling = ppg2.jobs._mark_function_wrapped(
            adapt_calling, calculating_function
        )
        job = ppg2.SharedMultiFileGeneratingJob(
            output_dir_prefix,
            output_files,
            adapt_calling,
            resources=resources,
            remove_unused=remove_unused,
            empty_ok=True,  # that's the default though
            # remove_build_dir_on_error=False, # if you ever need to debug the output
        )
        for input_file in input_files:
            job.depends_on(ppg.FileInvariant(input_file))
        for name, func in further_function_deps.items():
            job.depends_on(
                ppg.FunctionInvariant(str(output_dir_prefix) + "/" + name, func)
            )
        job.depends_on(ppg.ParameterInvariant(output_dir_prefix, version))
        job.version = version
        return job

    def _prebuild_ppg2_no_ppg(
        self,
        name,
        version,
        input_files,
        output_files,
        calculating_function,
        further_function_deps,
        resources,
        remove_unused,
        output_dir_prefix,
    ):
        import pypipegraph2 as ppg2
        import json

        history_filename = (
            Path(".ppg/history/") / ppg2.SharedMultiFileGeneratingJob.log_filename
        )
        if not history_filename.exists():
            raise ValueError(
                f"Could not find ppg history SharedMultiFileGeneratingJob log filename {history_filename}. Run a ppg with these objects before using them without ppg"
            )
        known = json.loads(history_filename.read_text())
        if not str(output_dir_prefix) in known:
            raise ValueError(
                f"Could not find {output_dir_prefix} in ppg history SharedMultiFileGeneratingJob log filename. Run a ppg with these objects before using them without pppg"
            )
        real_prefix = output_dir_prefix / "by_input" / known[str(output_dir_prefix)]
        if not real_prefix.exists():
            raise ValueError(
                "We could find a buld to use from the SharedMultiFileGeneratingJob log file - but it wasn't there"
            )
        for of in output_files:
            if not (real_prefix / of).exists():
                raise ValueError(
                    f"{real_prefix / of} was missing, but the folder was there. The most likely cause is that this has changed since the last ppg run and you need to run a ppg with these objects before using them outside of a ppg"
                )
        return DummyJob(real_prefix, [real_prefix / fn for fn in output_files])

    def _prebuild_ppg1(
        self,
        name,
        version,
        input_files,
        output_files,
        calculating_function,
        minimum_acceptable_version,
        maximum_acceptable_version,
        further_function_deps={},
    ):
        if minimum_acceptable_version is None:
            minimum_acceptable_version = version

        available_versions = self._find_versions(name)
        if version in available_versions:
            output_path = available_versions[version]
        else:
            # these are within minimum..maximum_acceptable_version
            acceptable_versions = sort_versions(
                [
                    (v, p)
                    for v, p in available_versions.items()
                    if (
                        (Version(v) >= minimum_acceptable_version)
                        and (
                            maximum_acceptable_version is None
                            or (Version(v) < maximum_acceptable_version)
                        )
                    )
                ]
            )
            ok_versions = []

            (
                new_source,
                new_funchash,
                new_closure,
            ) = ppg.FunctionInvariant._hash_function(calculating_function)

            for v, p in acceptable_versions:
                func_md5sum_path = p / "mbf.func.md5sum"
                func_md5sum_path2 = p / "mbf.func.md5sum2"
                try:
                    func_md5sum = json.loads(func_md5sum_path2.read_text())
                except OSError:
                    func_md5sum = func_md5sum_path.read_text()
                ok = False
                try:
                    ppg.FunctionInvariant._compare_new_and_old(
                        new_source, new_funchash, new_closure, func_md5sum
                    )
                    ok = False
                except ppg.NothingChanged:
                    ok = True
                if ok:
                    ok_versions.append((v, p))

            if ok_versions:
                version, output_path = ok_versions[-1]
            else:  # no version that is within the acceptable range and had the same build function
                output_path = self.prebuilt_path / self.hostname / name / version

        if ppg.inside_ppg():
            job = PrebuildJob(output_files, calculating_function, output_path)
            job.depends_on(_PrebuildFileInvariantsExploding(output_path, input_files))
            job.version = version
            return job
        else:
            for of in output_files:
                if not (output_path / of).exists():
                    raise ValueError(
                        "%s was missing and prebuild used outside of ppg - can't build it"
                        % (output_path / of).absolute()
                    )
            return DummyJob(
                output_path,
                PrebuildJob._normalize_output_files(output_files, output_path),
            )


_global_manager = None


def change_global_manager(new_manager):
    global _global_manager
    _global_manager = new_manager


def get_global_manager():
    return _global_manager


@contextlib.contextmanager
def with_global_manager(new_manager):
    """A context manager to be able to preplace the global manager for a block of code (tests...)"""
    old_manager = _global_manager
    change_global_manager(new_manager)
    try:
        yield
    finally:
        change_global_manager(old_manager)
