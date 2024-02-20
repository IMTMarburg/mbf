from pathlib import Path
import os
import time
import subprocess
from abc import ABC, abstractmethod
import pypipegraph as ppg
from .util import binary_exists


_external_algorithm_singletons = {}


def _hash_folder(path, hash):
    for fn in sorted(path.glob("*")):
        if fn.is_dir():
            _hash_folder(fn, hash)
        else:
            hash.update(fn.read_bytes())


def get_nix_store_path_from_binary(binary_name):
    bin = (Path("/bin") / binary_name).resolve()
    if not bin.exists():
        bin = Path(
            subprocess.check_output(["which", binary_name]).decode("utf-8")
        ).resolve()
    parts = bin.parts
    if parts[1] != "nix":
        raise ValueError(
            f"Something wrong with finding the nix name to {binary_name} from {bin}, {parts}"
        )
    nix_name = parts[3]
    return Path("/nix/store") / nix_name


def hash_folder(folder):
    """find the nix folder from that binary name,
    hash the complete folder"""
    import hashlib

    hash = hashlib.sha256()
    _hash_folder(folder, hash)
    hash = hash.hexdigest()
    return hash


class ExternalAlgorithm(ABC):
    """Together with an ExternalAlgorithmStore (or the global one),
    ExternalAlgorithm encapsulates a callable algorithm such as a high throughput aligner.

    They expect that there is at least a binary around in /bin
    """

    def __new__(cls):
        """For a given ExternalAlgorithm (by classname) only one object."""
        if not (cls.__name__,) in _external_algorithm_singletons:
            _external_algorithm_singletons[cls.__name__] = object.__new__(cls)
        return _external_algorithm_singletons[cls.__name__]

    def __init__(self):
        self.warn_if_binary_is_missing()
        self.version = self.get_version_cached()
        self.buffer_output = True

    def warn_if_binary_is_missing(self):
        if not binary_exists(self.primary_binary):
            import warnings

            warning = (
                f"Binary {self.primary_binary} not found.\n"
                "Try additing it via anysnake2.\n"
                "See https://github.com/IMTMarburg/flakes\n"
            )
            if hasattr(self, "flake_name"):
                warning += "(flake name is {self.flake_name})"
            warnings.warn(warning)

    def get_version(self):
        """
        This is here so you can overrite it e.g.
        with a call to /bin/algorithm --version

        The default turns the bin into a nix store path,
        and hashes that store path.

        The downstream ppg2 jobs rely on this
        """
        return hash_folder(get_nix_store_path_from_binary(self.primary_binary))

    def get_version_cached(self):
        """This makes sure you don't need to evaluate the version multiple times
        if it's the same nix store path
        """
        cache_path = Path("cache/nix-output-hashes")
        cache_path.mkdir(exist_ok=True, parents=True)
        nix_path = get_nix_store_path_from_binary(self.primary_binary)
        cache_fn = cache_path / nix_path.name
        if not cache_fn.exists():
            v = self.get_version()
            cache_fn.write_text(v)
            return v
        else:
            return cache_fn.read_text()

    @property
    @abstractmethod
    def name(self):
        pass  # pragma: no cover

    @property
    @abstractmethod
    def primary_binary(self):
        pass  # pragma: no cover

    @abstractmethod
    def build_cmd(self, output_directory, ncores, arguments):
        pass  # pragma: no cover

    @property
    def multi_core(self):
        return False

    def run(
        self,
        output_directory,
        arguments=None,
        cwd=None,
        call_afterwards=None,
        additional_files_created=None,
    ):
        """Return a job that runs the algorithm and puts the
        results in output_directory.
        Note that assigning different ouput_directories to different
        versions is your problem.
        """
        output_directory = Path(output_directory)
        output_directory.mkdir(parents=True, exist_ok=True)
        sentinel = output_directory / "sentinel.txt"
        filenames = [
            sentinel,
            output_directory / "stdout.txt",
            output_directory / "stderr.txt",
            output_directory / "cmd.txt",
        ]
        if additional_files_created:
            if isinstance(additional_files_created, (str, Path)):
                additional_files_created = [additional_files_created]
            filenames.extend(additional_files_created)

        job = ppg.MultiFileGeneratingJob(
            filenames,
            self.get_run_func(
                output_directory, arguments, cwd=cwd, call_afterwards=call_afterwards
            ),
            empty_ok=[output_directory / "stdout.txt", output_directory / "stderr.txt"],
        )
        job.depends_on(
            ppg.FunctionInvariant(job.job_id + "_call_afterwards", call_afterwards),
        )
        job.ignore_code_changes()
        job.depends_on(
            ppg.FunctionInvariant(
                job.job_id + "_build_cmd_func", self.__class__.build_cmd
            )
        )
        if self.multi_core:
            job.cores_needed = -1
            if hasattr(ppg, "is_ppg2"):
                import pypipegraph2 as ppg2

                assert job.resources == ppg2.Resources.AllCores

        return job

    def get_run_func(self, output_directory, arguments, cwd=None, call_afterwards=None):
        """Return the function that get_run_job passes to the job"""

        def do_run():
            sentinel = output_directory / "sentinel.txt"
            stdout = output_directory / "stdout.txt"
            stderr = output_directory / "stderr.txt"
            cmd_out = output_directory / "cmd.txt"

            op_stdout = open(stdout, "wb", buffering = -1 if self.buffer_output else 0)
            op_stderr = open(stderr, "wb", buffering = -1 if self.buffer_output else 0)
            cmd = [
                str(x)
                for x in self.build_cmd(
                    output_directory,
                    ppg.util.global_pipegraph.rc.cores_available
                    if self.multi_core
                    else 1,
                    arguments() if callable(arguments) else arguments,
                )
            ]
            cmd_out.write_text(" ".join(cmd))
            start_time = time.time()
            print(" ".join(cmd))
            env = os.environ.copy()
            if (
                "LD_LIBRARY_PATH" in env
            ):  # rpy2 likes to sneak this in, breaking e.g. STAR
                del env["LD_LIBRARY_PATH"]
            p = subprocess.Popen(
                cmd, stdout=op_stdout, stderr=op_stderr, cwd=cwd, env=env
            )
            p.communicate()
            op_stdout.close()
            op_stderr.close()
            ok = self.check_success(
                p.returncode, stdout.read_bytes(), stderr.read_bytes()
            )
            if ok is True:
                runtime = time.time() - start_time
                sentinel.write_text(
                    f"run time: {runtime:.2f} seconds\nreturn code: {p.returncode}"
                )
                if call_afterwards is not None:
                    call_afterwards()
            else:
                raise ValueError(
                    f"{self.name} run failed. Error was: {ok}. Cmd was: {cmd}"
                )

        return do_run

    def check_success(self, return_code, stdout, stderr):
        if return_code == 0:
            return True
        else:
            return f"Return code != 0: {return_code}"
