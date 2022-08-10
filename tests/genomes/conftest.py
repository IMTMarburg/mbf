#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import sys

# import pytest
from pathlib import Path

import pypipegraph2 as ppg2  # noqa: F401


ppg2.replace_ppg1()  # noqa: F401


from pypipegraph.testing.fixtures import (  # noqa: F401
    new_pipegraph,  # noqa: F401
    pytest_runtest_makereport,  # noqa: F401
)  # noqa: F401

root = Path(__file__).parent.parent
sys.path.append(str(root / "src"))

from mbf_genomes.testing.fixtures import mock_download, shared_prebuild  # noqa: F401
