from .star import STAR, STARSolo
from .subread import Subread, Subjunc
from .bowtie import Bowtie
from .salmon import Salmon
from .bwa import BWA
from .bbmap import BBMap

__all__ = ["Bowtie", "Subread", "STAR", "STARSolo", "Salmon", "Subjunc", "BWA", "BBMap"]

try:
    from .bbmap import ExtendCigarBBMap

    __all__.append(ExtendCigarBBMap)
except ImportError:
    pass
