"""Tests for the Milestone 2 synthetic datacube workflow."""

import numpy as np

from molecular_conditions.datacube import load_fits_cube, summarize_cube_metadata
from molecular_conditions.demo_data import create_synthetic_line_cube
import pytest

from molecular_conditions.moments import make_moment0, make_moment1, make_moment2
from molecular_conditions.spectra import (
    estimate_rms_noise,
    extract_pixel_spectrum,
    extract_region_spectrum,
)


def test_load_fits_cube_missing_path_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_fits_cube(tmp_path / "missing_cube.fits")


def test_synthetic_fits_cube_creation(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")

    assert cube_path.exists()
    assert cube_path.stat().st_size > 0


def test_cube_loading_and_metadata_summary(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)
    summary = summarize_cube_metadata(cube)

    expected_keys = {
        "shape",
        "unit",
        "spectral_axis_unit",
        "spectral_axis_min",
        "spectral_axis_max",
        "wcs_summary",
        "beam",
        "has_valid_spectral_axis",
    }
    assert expected_keys.issubset(summary)
    assert summary["shape"] == (64, 32, 32)
    assert summary["unit"] == "K"
    assert summary["has_valid_spectral_axis"] is True
    assert summary["wcs_summary"] is not None
    assert summary["beam"] is not None


def test_pixel_spectrum_extraction_shape(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)

    spectral_axis, intensity = extract_pixel_spectrum(cube, x=16, y=15)

    assert spectral_axis.shape == (64,)
    assert intensity.shape == (64,)
    assert str(intensity.unit) == "K"


def test_pixel_spectrum_out_of_bounds_raises_index_error(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)

    with pytest.raises(IndexError):
        extract_pixel_spectrum(cube, x=32, y=15)


def test_region_spectrum_extraction_shape_and_unit(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)
    region = {"x_min": 12, "x_max": 20, "y_min": 11, "y_max": 19}

    spectral_axis, intensity = extract_region_spectrum(cube, region, statistic="mean")

    assert spectral_axis.shape == (64,)
    assert intensity.shape == (64,)
    assert str(intensity.unit) == "K"


def test_rms_noise_positive_and_finite(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)
    spectrum = extract_pixel_spectrum(cube, x=16, y=15)

    rms = estimate_rms_noise(spectrum, line_free_window=(-6300.0, -4300.0))

    assert np.isfinite(rms)
    assert rms > 0.0


def test_rms_noise_can_return_quantity(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)
    spectrum = extract_pixel_spectrum(cube, x=16, y=15)

    rms = estimate_rms_noise(
        spectrum,
        line_free_window=(-6300.0, -4300.0),
        return_quantity=True,
    )

    assert np.isfinite(rms.value)
    assert rms.value > 0.0
    assert str(rms.unit) == "K"


def test_rms_noise_non_overlapping_window_raises_value_error(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)
    spectrum = extract_pixel_spectrum(cube, x=16, y=15)

    with pytest.raises(ValueError, match="does not overlap"):
        estimate_rms_noise(spectrum, line_free_window=(100_000.0, 101_000.0))


def test_moment_maps_have_expected_spatial_dimensions(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)

    moment0 = make_moment0(cube, spectral_window=(-3500.0, 3500.0))
    moment1 = make_moment1(cube, spectral_window=(-3500.0, 3500.0))
    moment2 = make_moment2(cube, spectral_window=(-3500.0, 3500.0))

    assert moment0.shape == (32, 32)
    assert moment1.shape == (32, 32)
    assert moment2.shape == (32, 32)
    assert np.isfinite(moment0.value).any()
    assert np.isfinite(moment1.value).any()
    assert np.isfinite(moment2.value).any()


def test_moment0_accepts_reversed_spectral_window(tmp_path):
    cube_path = create_synthetic_line_cube(tmp_path / "synthetic_line_cube.fits")
    cube = load_fits_cube(cube_path)

    moment0 = make_moment0(cube, spectral_window=(3500.0, -3500.0))

    assert moment0.shape == (32, 32)
    assert np.isfinite(moment0.value).any()
