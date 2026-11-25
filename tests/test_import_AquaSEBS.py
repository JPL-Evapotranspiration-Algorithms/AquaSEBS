"""Smoke tests for top-level AquaSEBS package imports and metadata."""


def test_import_AquaSEBS():
    import AquaSEBS

    assert hasattr(AquaSEBS, "AquaSEBS")
    assert hasattr(AquaSEBS, "water_heat_flux")
    assert hasattr(AquaSEBS, "__version__")
    assert hasattr(AquaSEBS, "__author__")
