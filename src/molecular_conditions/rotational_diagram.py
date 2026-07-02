"""LTE rotational diagram analysis.

The implementation here is intentionally compact and explicit. It assumes
optically thin LTE emission, a single excitation temperature, beam filling
factor near one, and line measurements from one molecular species.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy import constants as const


REQUIRED_LINE_COLUMNS = {
    "species",
    "transition",
    "rest_frequency_ghz",
    "upper_energy_k",
    "einstein_a_s",
    "degeneracy_upper",
    "integrated_intensity_K_km_s",
}

OPTIONAL_ERROR_COLUMN = "integrated_intensity_error_K_km_s"

DEFAULT_ASSUMPTIONS = (
    "Optically thin LTE emission; single rotational temperature; "
    "beam filling factor approximately one; same species and emitting region; "
    "column density is reported as N_tot / Q(T_rot), not absolute N_tot."
)


@dataclass(frozen=True)
class RotationalDiagramFit:
    """LTE rotational diagram fit result."""

    excitation_temperature_k: float
    intercept: float
    slope: float
    column_density_over_partition_function_cm2: float
    excitation_temperature_error_k: float | None = None
    intercept_error: float | None = None
    slope_error: float | None = None
    column_density_over_partition_function_error_cm2: float | None = None
    chi_square: float | None = None
    reduced_chi_square: float | None = None
    assumptions: str | list[str] = field(default=DEFAULT_ASSUMPTIONS)

    def as_dict(self) -> dict[str, float | str | list[str] | None]:
        """Return fit results in a table/export-friendly dictionary."""
        return {
            "excitation_temperature_k": self.excitation_temperature_k,
            "excitation_temperature_error_k": self.excitation_temperature_error_k,
            "intercept": self.intercept,
            "intercept_error": self.intercept_error,
            "slope": self.slope,
            "slope_error": self.slope_error,
            "column_density_over_partition_function_cm2": (
                self.column_density_over_partition_function_cm2
            ),
            "column_density_over_partition_function_error_cm2": (
                self.column_density_over_partition_function_error_cm2
            ),
            "chi_square": self.chi_square,
            "reduced_chi_square": self.reduced_chi_square,
            "assumptions": self.assumptions,
        }


def _load_line_measurements(line_measurements: str | Path | pd.DataFrame) -> pd.DataFrame:
    """Load line measurements from a path or copy an existing DataFrame."""
    if isinstance(line_measurements, (str, Path)):
        path = Path(line_measurements)
        if not path.exists():
            raise FileNotFoundError(f"Line-measurement table does not exist: {path}")
        return pd.read_csv(path, comment="#")
    if isinstance(line_measurements, pd.DataFrame):
        return line_measurements.copy()
    raise TypeError("line_measurements must be a CSV path or pandas DataFrame.")


def _validate_line_measurements(table: pd.DataFrame) -> pd.DataFrame:
    """Validate required columns, same-species constraint, and positive values."""
    missing_columns = REQUIRED_LINE_COLUMNS - set(table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Line-measurement table is missing required columns: {missing}")

    table = table.copy()
    numeric_columns = [
        "rest_frequency_ghz",
        "upper_energy_k",
        "einstein_a_s",
        "degeneracy_upper",
        "integrated_intensity_K_km_s",
    ]
    if OPTIONAL_ERROR_COLUMN in table.columns:
        numeric_columns.append(OPTIONAL_ERROR_COLUMN)
    for column in numeric_columns:
        table[column] = pd.to_numeric(table[column], errors="raise")

    species = table["species"].astype(str).unique()
    if len(species) != 1:
        raise ValueError(
            "Rotational diagrams must use transitions from one species only; "
            f"found {list(species)}."
        )

    positive_columns = [
        "rest_frequency_ghz",
        "einstein_a_s",
        "degeneracy_upper",
        "integrated_intensity_K_km_s",
    ]
    if OPTIONAL_ERROR_COLUMN in table.columns:
        positive_columns.append(OPTIONAL_ERROR_COLUMN)
    for column in positive_columns:
        if (table[column] <= 0).any():
            raise ValueError(f"{column} must contain only positive values.")

    return table.sort_values("upper_energy_k").reset_index(drop=True)


def build_population_diagram(line_measurements: str | Path | pd.DataFrame) -> pd.DataFrame:
    """Convert integrated intensities into rotational-diagram points.

    For optically thin emission in brightness-temperature units, the upper
    state column density is

    ``N_u = 8*pi*k*nu^2 / (h*c^3*A_ul) * integral(T_B dv)``.

    Frequencies are read from GHz, integrated intensities from K km/s, and the
    returned ``N_u`` values are in ``cm^-2``. The logarithmic diagram quantity
    is ``ln(N_u / g_u)``.
    """
    table = _validate_line_measurements(_load_line_measurements(line_measurements))

    frequency_hz = table["rest_frequency_ghz"].to_numpy(dtype=float) * 1.0e9
    integrated_intensity_k_m_s = (
        table["integrated_intensity_K_km_s"].to_numpy(dtype=float) * 1000.0
    )
    einstein_a_s = table["einstein_a_s"].to_numpy(dtype=float)
    degeneracy_upper = table["degeneracy_upper"].to_numpy(dtype=float)

    factor_m2 = (
        8.0
        * np.pi
        * const.k_B.value
        * frequency_hz**2
        / (const.h.value * const.c.value**3 * einstein_a_s)
    )
    upper_column_density_cm2 = factor_m2 * integrated_intensity_k_m_s / 1.0e4
    upper_column_density_over_degeneracy_cm2 = (
        upper_column_density_cm2 / degeneracy_upper
    )
    ln_nu_over_gu = np.log(upper_column_density_over_degeneracy_cm2)

    population = table.copy()
    population["upper_column_density_cm2"] = upper_column_density_cm2
    population["upper_column_density_over_degeneracy_cm2"] = (
        upper_column_density_over_degeneracy_cm2
    )
    population["ln_Nu_over_gu"] = ln_nu_over_gu

    if OPTIONAL_ERROR_COLUMN in population.columns:
        intensity = population["integrated_intensity_K_km_s"].to_numpy(dtype=float)
        intensity_error = population[OPTIONAL_ERROR_COLUMN].to_numpy(dtype=float)
        fractional_error = intensity_error / intensity
        population["upper_column_density_error_cm2"] = (
            upper_column_density_cm2 * fractional_error
        )
        population["ln_Nu_over_gu_error"] = fractional_error

    return population


def _linear_fit_with_covariance(
    x: np.ndarray,
    y: np.ndarray,
    y_error: np.ndarray | None,
) -> tuple[float, float, np.ndarray, float | None, float | None]:
    """Fit ``y = intercept + slope*x`` and return parameters and covariance."""
    design = np.column_stack([np.ones_like(x), x])
    dof = max(0, x.size - 2)

    if y_error is not None:
        weights = 1.0 / y_error**2
        weighted_design = design * weights[:, None]
        normal_matrix = design.T @ weighted_design
        covariance = np.linalg.inv(normal_matrix)
        beta = covariance @ (design.T @ (weights * y))
        residual = y - design @ beta
        chi_square = float(np.sum((residual / y_error) ** 2))
        reduced_chi_square = chi_square / dof if dof > 0 else None
    else:
        beta, residual_sum, _, _ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design @ beta
        chi_square = float(np.sum(residual**2))
        reduced_chi_square = chi_square / dof if dof > 0 else None
        covariance = np.linalg.inv(design.T @ design)
        if dof > 0:
            covariance = covariance * reduced_chi_square
        if residual_sum.size:
            chi_square = float(residual_sum[0])

    intercept, slope = [float(value) for value in beta]
    return intercept, slope, covariance, chi_square, reduced_chi_square


def fit_rotational_temperature(population_table: pd.DataFrame) -> RotationalDiagramFit:
    """Fit LTE rotational temperature from a population diagram.

    The fitted relation is ``ln(N_u / g_u) = intercept + slope * E_u`` where
    ``E_u`` is in Kelvin. Therefore ``T_rot = -1 / slope``. The intercept is
    ``ln(N_tot / Q(T_rot))``; without a partition function, this function
    reports ``N_tot / Q(T_rot)`` rather than an absolute total column density.
    """
    if not isinstance(population_table, pd.DataFrame):
        raise TypeError("population_table must be a pandas DataFrame.")

    required_columns = {"upper_energy_k", "ln_Nu_over_gu"}
    missing_columns = required_columns - set(population_table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Population table is missing required columns: {missing}")

    x = population_table["upper_energy_k"].to_numpy(dtype=float)
    y = population_table["ln_Nu_over_gu"].to_numpy(dtype=float)
    if x.size < 3:
        raise ValueError("At least three transitions are required for a robust fit.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Population diagram columns must contain finite values.")
    if np.nanmax(x) == np.nanmin(x):
        raise ValueError("upper_energy_k must span a nonzero range.")

    y_error = None
    if "ln_Nu_over_gu_error" in population_table.columns:
        candidate = population_table["ln_Nu_over_gu_error"].to_numpy(dtype=float)
        if not np.all(np.isfinite(candidate)):
            raise ValueError(
                "ln_Nu_over_gu_error must contain only finite values when provided."
            )
        if not np.all(candidate > 0):
            raise ValueError(
                "ln_Nu_over_gu_error must contain only positive values when provided."
            )
        y_error = candidate

    intercept, slope, covariance, chi_square, reduced_chi_square = (
        _linear_fit_with_covariance(x, y, y_error)
    )
    if slope >= 0:
        raise ValueError(
            "Fitted rotational-diagram slope is non-negative; cannot infer a "
            "positive rotational temperature."
        )

    excitation_temperature_k = -1.0 / slope
    column_density_over_partition_function_cm2 = float(np.exp(intercept))

    intercept_error = None
    slope_error = None
    excitation_temperature_error_k = None
    column_density_over_partition_function_error_cm2 = None
    diagonal = np.diag(covariance)
    if np.all(np.isfinite(diagonal)) and np.all(diagonal >= 0):
        intercept_error = float(np.sqrt(diagonal[0]))
        slope_error = float(np.sqrt(diagonal[1]))
        excitation_temperature_error_k = float(slope_error / slope**2)
        column_density_over_partition_function_error_cm2 = float(
            column_density_over_partition_function_cm2 * intercept_error
        )

    return RotationalDiagramFit(
        excitation_temperature_k=float(excitation_temperature_k),
        intercept=intercept,
        slope=slope,
        column_density_over_partition_function_cm2=(
            column_density_over_partition_function_cm2
        ),
        excitation_temperature_error_k=excitation_temperature_error_k,
        intercept_error=intercept_error,
        slope_error=slope_error,
        column_density_over_partition_function_error_cm2=(
            column_density_over_partition_function_error_cm2
        ),
        chi_square=chi_square,
        reduced_chi_square=reduced_chi_square,
    )
