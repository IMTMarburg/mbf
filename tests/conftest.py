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


# from pypipegraph.testing.fixtures import new_pipegraph, pytest_runtest_makereport  # noqa:F401

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))


def pytest_generate_tests(metafunc):
    if "both_ppg_and_no_ppg" in metafunc.fixturenames:
        metafunc.parametrize("both_ppg_and_no_ppg", [True, False], indirect=True)
