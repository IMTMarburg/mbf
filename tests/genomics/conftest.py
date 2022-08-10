#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import sys
import pytest  # noqa:F401
from pathlib import Path
import pypipegraph2  # noqa:F401

pypipegraph2.replace_ppg1()

from pypipegraph.testing.fixtures import (  # noqa:F401
    new_pipegraph,
    both_ppg_and_no_ppg,
    no_pipegraph,
    pytest_runtest_makereport,
)  # noqa:F401
from pypipegraph2.testing.fixtures import job_trace_log  # noqa:F401
from mbf_qualitycontrol.testing.fixtures import (  # noqa:F401
    new_pipegraph_no_qc,
    both_ppg_and_no_ppg_no_qc,
)
from mbf_genomics.testing.fixtures import clear_annotators  # noqa:F401

# from mbf_externals.testing.fixtures import local_store, global_store  # noqa:F401
root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))

from plotnine.tests.conftest import (  # noqa:F401
    _setup,
    _teardown,  # noqa:F401
    pytest_assertrepr_compare,
)

_setup()


def pytest_generate_tests(metafunc):
    if "both_ppg_and_no_ppg" in metafunc.fixturenames:
        metafunc.parametrize("both_ppg_and_no_ppg", [True, False], indirect=True)
