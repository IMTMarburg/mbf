from pathlib import Path
from .base import GenomeBase, HardCodedGenome
from .ensembl import EnsemblGenome
from .filebased import FileBasedGenome, InteractiveFileBasedGenome
from .igenomes import IGenome

data_path = Path(__file__).parent / "data"
__version__ = "0.4"


def Homo_sapiens(rev):
    return EnsemblGenome("Homo_sapiens", rev)


def Mus_musculus(rev):
    return EnsemblGenome("Mus_musculus", rev)


__all__ = [
    "GenomeBase",
    "EnsemblGenome",
    "FileBasedGenome",
    "InteractiveFileBasedGenome",
    "data_path",
    "HardCodedGenome",
    "Homo_sapiens",
    "Mus_musculus",
]
