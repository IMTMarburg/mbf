#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for mbf_align.

    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import pathlib

import pypipegraph2  # noqa:F401

pypipegraph2.replace_ppg1()  # noqa:F401


from pypipegraph.testing.fixtures import (  # noqa:F401
    new_pipegraph,
    pytest_runtest_makereport,
)

from mbf_qualitycontrol.testing.fixtures import new_pipegraph_no_qc  # noqa:F401

root = pathlib.Path(__file__).parent.parent


# sys.path.append(str(root / "src"))
# print("the path is", sys.path)
# subprocess.check_call(["python3", "setup.py", "build_ext", "-i"], cwd=root)
