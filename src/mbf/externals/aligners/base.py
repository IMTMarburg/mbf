from ..externals import ExternalAlgorithm
from pathlib import Path
from abc import abstractmethod
import pypipegraph as ppg


class Aligner(ExternalAlgorithm):
    @abstractmethod
    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_basename,
        output_bam_filename,
        parameters,
    ):
        pass  # pragma: no cover

    @abstractmethod
    def build_index_func(self, fasta_files, gtf_input_filename, output_prefix):
        pass  # pragma: no cover

    @abstractmethod
    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        pass  # pragma: no cover

    def build_cmd(self, output_dir, ncores, arguments):
        if (
            not isinstance(arguments, list)
            or len(arguments) < 2
            or arguments[0] != "FROM_ALIGNER"
        ):
            raise ValueError(
                "Please call one of the following functions instead: .align_job, .build_index_job"
                + str(arguments)
            )
        return self._aligner_build_cmd(output_dir, ncores, arguments[1:])

    def build_index_job(
        self, fasta_files, gtf_input_filename, output_fileprefix
    ):  # todo: is this used *outside* of test cases?
        output_directory = Path(output_fileprefix)
        output_directory.mkdir(parents=True, exist_ok=True)
        sentinel = output_directory / "sentinel.txt"
        job = ppg.FileGeneratingJob(
            sentinel,
            self.build_index_func(fasta_files, gtf_input_filename, output_directory),
        )
        if self.multi_core:
            job.cores_needed = -1
        job.index_path = output_fileprefix
        return job

    def build_index(self, fasta_files, gtf_input_filename, output_fileprefix):
        output_fileprefix = Path(output_fileprefix)
        output_fileprefix.mkdir(parents=True, exist_ok=True)
        func = self.build_index_func(fasta_files, gtf_input_filename, output_fileprefix)
        print(func)
        print(func.__code__)
        func()

    def get_index_version_range(self):  # pragma: no cover
        return None, None

    def get_index_filenames(self):
        """A list of filenames that build_index will create."""
        raise NotImplementedError()
