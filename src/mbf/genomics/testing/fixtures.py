import pytest


@pytest.fixture
def clear_annotators(request):
    """Clear the annotator singleton instance cache
    which is only used if no ppg is in play"""
    import mbf.genomics.annotator

    mbf.genomics.annotator.annotator_singletons.clear()
    mbf.genomics.annotator.annotator_singletons["lookup"] = []
