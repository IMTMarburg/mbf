import pytest


@pytest.fixture
def clear_annotators(request):
    """Clear the annotator singleton instance cache
    which is only used if no ppg is in play"""
    import mbf_genomics.annotator

    mbf_genomics.annotator.annotator_singletons.clear()
    mbf_genomics.annotator.annotator_singletons["lookup"] = []
