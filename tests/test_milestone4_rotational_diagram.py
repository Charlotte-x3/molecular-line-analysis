"""Tests for Milestone 4 LTE rotational diagram analysis."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from molecular_conditions.rotational_diagram import (
    build_population_diagram,
    fit_rotational_temperature,
)


DEMO_ROTATIONAL_LINES = Path("data/demo/demo_rotational_lines.csv")


def test_build_population_diagram_returns_expected_columns():
    population = build_population_diagram(DEMO_ROTATIONAL_LINES)

    expected_columns = {
        "species",
        "transition",
        "rest_frequency_ghz",
        "upper_energy_k",
        "einstein_a_s",
        "degeneracy_upper",
        "integrated_intensity_K_km_s",
        "upper_column_density_cm2",
        "upper_column_density_over_degeneracy_cm2",
        "ln_Nu_over_gu",
        "ln_Nu_over_gu_error",
    }
    assert expected_columns.issubset(population.columns)
    assert len(population) >= 5


def test_population_diagram_values_are_finite_and_positive():
    population = build_population_diagram(DEMO_ROTATIONAL_LINES)

    assert np.isfinite(population["upper_column_density_cm2"]).all()
    assert np.isfinite(population["upper_column_density_over_degeneracy_cm2"]).all()
    assert np.isfinite(population["ln_Nu_over_gu"]).all()
    assert (population["upper_column_density_cm2"] > 0).all()
    assert (population["upper_column_density_over_degeneracy_cm2"] > 0).all()


def test_fit_rotational_temperature_recovers_synthetic_temperature():
    population = build_population_diagram(DEMO_ROTATIONAL_LINES)
    fit = fit_rotational_temperature(population)

    assert fit.excitation_temperature_k == pytest.approx(38.0, rel=0.15)
    assert fit.column_density_over_partition_function_cm2 > 0.0
    assert fit.slope < 0.0
    assert fit.as_dict()["excitation_temperature_k"] == fit.excitation_temperature_k


def test_fit_rotational_temperature_without_uncertainties():
    population = build_population_diagram(DEMO_ROTATIONAL_LINES)
    population = population.drop(columns=["ln_Nu_over_gu_error"])

    fit = fit_rotational_temperature(population)

    assert fit.excitation_temperature_k == pytest.approx(38.0, rel=0.15)
    assert fit.reduced_chi_square is not None


def test_invalid_population_uncertainties_raise_value_error():
    population = build_population_diagram(DEMO_ROTATIONAL_LINES)
    population.loc[0, "ln_Nu_over_gu_error"] = 0.0

    with pytest.raises(ValueError, match="positive"):
        fit_rotational_temperature(population)


def test_mixed_species_raise_clear_value_error():
    measurements = pd.read_csv(DEMO_ROTATIONAL_LINES, comment="#")
    measurements.loc[0, "species"] = "OtherMol"

    with pytest.raises(ValueError, match="one species"):
        build_population_diagram(measurements)


def test_missing_required_columns_raise_clear_value_error():
    measurements = pd.read_csv(DEMO_ROTATIONAL_LINES, comment="#")
    measurements = measurements.drop(columns=["einstein_a_s"])

    with pytest.raises(ValueError, match="missing required columns"):
        build_population_diagram(measurements)
