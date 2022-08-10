import pypipegraph as ppg

__version__ = "0.1"


def register_qc(job):
    for k in job:
        if not isinstance(job, ppg.Job):
            raise TypeError("register_qc takes only job objects")
    job._mbf_qc = True
    if hasattr(ppg.util.global_pipegraph, "_qc_keep_function") and (
        getattr(ppg.util.global_pipegraph, "_qc_keep_function") is False
        or not getattr(ppg.util.global_pipegraph, "_qc_keep_function")(job)
    ):
        job.prune()
    for attr in ["lfg", "cache_job", "table_job"]:
        j = getattr(job, attr, None)
        if j is not None:
            if isinstance(j, ppg.Job):
                register_qc(j)
    return job


def qc_disabled():
    if not ppg.inside_ppg():
        return True
    return getattr(ppg.util.global_pipegraph, "_qc_keep_function", True) is False


def disable_qc():
    """Disable qc.
    That means new jobs that are generated are automatically pruned
    (but may be revived by calling prune_qc with an appropriate keep function),
    but code that depends on qc_disabled() does not even generate the jobs"""

    ppg.util.global_pipegraph._qc_keep_function = False


def prune_qc(keep=False):
    """Prune all qc jobs but those where keep returns True.
    Also applies to further qc jobs!"""
    ppg.util.global_pipegraph._qc_keep_function = keep
    for job in get_qc_jobs():
        if keep is not False and keep(job):
            job.unprune()
        else:
            job.prune()


def get_qc_jobs():
    """Get all qc jobs defined so far"""
    for job in ppg.util.global_pipegraph.jobs.values():
        if hasattr(job, "_mbf_qc"):
            yield job


class QCCollectingJob(ppg.FileGeneratingJob):
    def __init__(self, job_id, callback, depends_on_function=True):
        # if ppg.job.was_inited_before(self, QCCollectingJob):
        if hasattr(self, "inner_callback"):
            return
        self.inner_callback = callback

        def cb(output_filename):
            callback(output_filename, self.objects)

        self.objects = []
        if hasattr(ppg, "is_ppg2"):
            import pypipegraph2 as ppg2

            cb = ppg2.jobs._mark_function_wrapped(cb, self.inner_callback)
            super().__init__(job_id, cb)
            if not depends_on_function:
                self._ignore_code_changes()
            else:
                self.do_ignore_code_changes = False
                self.depends_on(self._create_parameter_dependency)
                self._handle_function_dependency(self.generating_function)
        else:
            super().__init__(job_id, cb)

    def add(self, obj):
        self.objects.append(obj)
        return self

    def inject_auto_invariants(self):
        if not self.do_ignore_code_changes:
            self.depends_on(
                ppg.FunctionInvariant(self.job_id + "_func", self.inner_callback)
            )
            self._create_parameter_dependency()

    def _create_parameter_dependency(self):
        names = []
        for obj in self.objects:
            if hasattr(obj, "name"):
                names.append(obj.name)
            elif hasattr(obj, "columns"):
                names.append(obj.columns[0])
            else:
                print(type(obj))
                raise ValueError(dir(obj))
        self.depends_on(
            ppg.ParameterInvariant(
                self.job_id,
                tuple(sorted(names)),
            )
        )
