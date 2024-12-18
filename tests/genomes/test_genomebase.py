import pytest
import pandas as pd

# from pathlib import Path
# from mbf.externals import PrebuildManager
# from mbf.genomes import EnsemblGenome
# import pypipegraph as ppg


@pytest.mark.usefixtures("new_pipegraph")
class TestBase:
    pass


def test_msgpack_unpacking_class_wrong_property_name():
    from mbf.genomes.base import msgpack_unpacking_class, MsgPackProperty

    with pytest.raises(NotImplementedError):

        @msgpack_unpacking_class
        class Shu:
            def _pepare_dfnolower():
                return pd.DataFrame()

            dfnolower = MsgPackProperty()
