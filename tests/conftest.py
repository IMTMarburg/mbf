#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import pytest
import sys
from pathlib import Path
import pypipegraph2 as ppg2  # noqa: F401


ppg2.replace_ppg1()

from pypipegraph.testing.fixtures import (  # noqa:F401
    new_pipegraph,
    no_pipegraph,
    both_ppg_and_no_ppg,
    pytest_runtest_makereport,
)
from mbf.qualitycontrol.testing.fixtures import (  # noqa: F401
    new_pipegraph_no_qc,
    both_ppg_and_no_ppg_no_qc,
)

from pypipegraph2.testing.fixtures import job_trace_log  # noqa:F401
from mbf.genomics.testing.fixtures import clear_annotators  # noqa:F401
from mbf.genomes.testing.fixtures import mock_download, shared_prebuild  # noqa: F40
import pytest
import os


# from pypipegraph.testing.fixtures import new_pipegraph, pytest_runtest_makereport  # noqa:F401

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))


def pytest_generate_tests(metafunc):
    if "both_ppg_and_no_ppg" in metafunc.fixturenames:
        metafunc.parametrize("both_ppg_and_no_ppg", [True, False], indirect=True)


@pytest.fixture()
def use_prebuild_genome(no_pipegraph):  # noqa: F811
    """Use a genome that has been previously build via ppg.
    This is a hacky hack, but needs must.

    Run a ppg with

    genome = mbf.genomes.Homo_sapiens(108)
    genome_mouse = mbf.genomes.Mus_musculus(108)
    compara = EnsemblCompara(108)
    genome.download()
    genome_mouse.download()
    compara.download()
    if it's not build on this machine.
    (
        if you're missing a version you'll need to merge the .ppg/history/SharedMultiFileGeneratingJobs.json to this folder
    )
    """
    if not "/run" in os.getcwd():
        raise ValueError(
            "Must be used in conjecture with no_pipegraph - path was ", str(__file__)
        )
    fn = Path(__file__).parent / "SharedMultiFileGeneratingJobs.json"
    input = fn.read_text()
    nice_hostname = os.environ.get("NICE_HOSTNAME")
    if nice_hostname:
        input = input.replace("/clara/", f"/{nice_hostname}/")
    Path(".ppg/history").mkdir(exist_ok=True, parents=True)
    (Path(".ppg/history") / "SharedMultiFileGeneratingJobs.json").write_text(input)
