import pytest
import pypipegraph as ppg


@pytest.fixture
def new_pipegraph_no_qc(new_pipegraph):
    ppg.util.global_pipegraph._qc_keep_function = False
    return new_pipegraph
    # this really does not work :(


@pytest.fixture
def both_ppg_and_no_ppg_no_qc(both_ppg_and_no_ppg):
    if ppg.util.global_pipegraph is not None:
        ppg.util.global_pipegraph._qc_keep_function = False
    return both_ppg_and_no_ppg
