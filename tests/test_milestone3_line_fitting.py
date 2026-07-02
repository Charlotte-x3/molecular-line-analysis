"""Tests for Milestone 3 transition matching and Gaussian fitting."""

from pathlib import Path

import numpy as np
import pytest
from astropy import units as u

from molecular_conditions.datacube import load_fits_cube
from molecular_conditions.demo_data import create_synthetic_line_cube
from molecular_conditions.fitting import (
    fit_single_gaussian,
    gaussian_model,
    integrated_intensity_from_gaussian,
)
from molecular_conditions.spectra import extract_pixel_spectrum
from molecular_conditions.transitions import load_transition_table, match_transitions


DEMO_CATALOG = Path("data/demo/transitions_demo.csv")


def test_load_transition_table_loads_demo_catalog():
    catalog = load_transition_table(DEMO_CATALOG)

    assert len(catalog) >= 5
    assert {
        "species",
        "quantum_numbers",
        "rest_frequency_ghz",
        "upper_energy_k",
        "einstein_a_s",
        "degeneracy_upper",
    }.issubset(catalog.columns)


def test_match_transitions_finds_co_10_near_115_271_ghz():
    catalog = load_transition_table(DEMO_CATALOG)
    matches = match_transitions(115.2712018, catalog, tolerance_mhz=5.0)

    assert len(matches) == 1
    assert matches.loc[0, "species"] == "CO"
    assert matches.loc[0, "quantum_numbers"] == "J=1-0"
    assert matches.loc[0, "abs_frequency_offset_mhz"] < 1.0e-6


def test_fit_single_gaussian_recovers_synthetic_parameters():
    rng = np.random.default_rng(10)
    x = np.linspace(-5.0, 5.0, 121)
    true_amplitude = 4.2
    true_centroid = 0.75
    true_sigma = 0.9
    true_baseline = 0.15
    y = gaussian_model(x, true_amplitude, true_centroid, true_sigma, true_baseline)
    y = y + rng.normal(0.0, 0.03, size=x.shape)

    fit = fit_single_gaussian(x, y)

    assert fit.amplitude == pytest.approx(true_amplitude, rel=0.05)
    assert fit.centroid == pytest.approx(true_centroid, abs=0.05)
    assert fit.sigma == pytest.approx(true_sigma, rel=0.08)
    assert fit.baseline == pytest.approx(true_baseline, abs=0.08)
    assert fit.integrated_intensity == pytest.approx(
        integrated_intensity_from_gaussian(fit.amplitude, fit.sigma)
    )


def test_fit_single_gaussian_emission_line_has_positive_width_and_integral():
    x = np.linspace(-4.0, 4.0, 80)
    y = gaussian_model(x, 2.0, 0.0, 0.7, 0.05)

    fit = fit_single_gaussian(x, y)

    assert fit.sigma > 0.0
    assert fit.integrated_intensity > 0.0
    assert fit.fwhm > 0.0


def test_fit_single_gaussian_records_quantity_units():
    x = np.linspace(-4.0, 4.0, 80) * (u.km / u.s)
    y = gaussian_model(x.value, 2.0, 0.0, 0.7, 0.05) * u.K

    fit = fit_single_gaussian(x, y)
    fit_dict = fit.as_dict()

    assert fit.x_unit == "km / s"
    assert fit.intensity_unit == "K"
    assert fit.integrated_intensity_unit == "K km / s"
    assert fit_dict["x_unit"] == "km / s"
    assert fit_dict["intensity_unit"] == "K"
    assert fit_dict["integrated_intensity_unit"] == "K km / s"
    assert fit_dict["fwhm"] == pytest.approx(fit.fwhm)


def test_fit_single_gaussian_invalid_inputs_raise_clear_errors():
    x = np.linspace(-1.0, 1.0, 10)
    y = np.ones(10)

    with pytest.raises(ValueError, match="nonzero variation"):
        fit_single_gaussian(x, y)

    with pytest.raises(ValueError, match="same shape"):
        fit_single_gaussian(x, y[:-1])

    y_with_nan = np.arange(10.0)
    y_with_nan[4] = np.nan
    with pytest.raises(ValueError, match="finite"):
        fit_single_gaussian(x, y_with_nan)


def test_gaussian_fitting_workflow_on_synthetic_cube_spectrum(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)
    spectral_axis, intensity = extract_pixel_spectrum(cube, x=16, y=15)

    fit = fit_single_gaussian(spectral_axis.to(u.km / u.s), intensity)

    assert fit.amplitude > 3.0
    assert abs(fit.centroid) < 0.5
    assert 0.4 < fit.sigma < 1.2
    assert fit.integrated_intensity > 0.0
    assert fit.x_unit == "km / s"
    assert fit.intensity_unit == "K"
    assert fit.integrated_intensity_unit == "K km / s"
