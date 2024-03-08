import os
from pathlib import Path
from .externals import (
    ExternalAlgorithm,
)
from .fastq import FASTQC
from .prebuild import PrebuildManager, change_global_manager, get_global_manager, with_global_manager
from . import aligners, peak_callers
from . import util


__version__ = "0.4"


def create_defaults():

    if "MBF_EXTERNAL_PREBUILD_PATH" in os.environ:
        hostname = os.environ["MBF_EXTERNAL_HOSTNAME"]
        if not (Path(os.environ["MBF_EXTERNAL_PREBUILD_PATH"]) / hostname).exists():
            raise ValueError(
                "%s did not exist - must be created manually"
                % (Path(os.environ["MBF_EXTERNAL_PREBUILD_PATH"]) / hostname)
            )
        prebuild_path = Path(os.environ["MBF_EXTERNAL_PREBUILD_PATH"])

    elif "VIRTUAL_ENV" in os.environ:
        import socket

        prebuild_path = (Path(".") / "prebuilt").absolute()
        prebuild_path.mkdir(exist_ok=True)
        hostname = socket.gethostname()
    else:
        # print("No defaults for mbf.externals possible")
        change_global_manager(None)
        return

    change_global_manager(PrebuildManager(prebuild_path, hostname))


create_defaults()


__all__ = [
    ExternalAlgorithm,
    FASTQC,
    PrebuildManager,
    aligners,
    peak_callers,
    util,
    __version__,
    create_defaults,
    change_global_manager,
    get_global_manager,
]
