"""Basic package import tests for the initial skeleton."""


def test_package_imports():
    import molecular_conditions

    assert molecular_conditions.__version__ == "0.1.0"


def test_placeholder_modules_import():
    from molecular_conditions import datacube
    from molecular_conditions import fitting
    from molecular_conditions import moments
    from molecular_conditions import radex_grid
    from molecular_conditions import rotational_diagram
    from molecular_conditions import spectra
    from molecular_conditions import transitions

    assert datacube is not None
    assert spectra is not None
    assert moments is not None
    assert fitting is not None
    assert transitions is not None
    assert rotational_diagram is not None
    assert radex_grid is not None
