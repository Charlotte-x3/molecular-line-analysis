"""Gaussian molecular-line fitting tools."""

from dataclasses import dataclass

import numpy as np
from astropy import units as u
from scipy.optimize import curve_fit


@dataclass(frozen=True)
class GaussianLineFit:
    """Summary of a fitted Gaussian spectral component."""

    amplitude: float
    centroid: float
    sigma: float
    integrated_intensity: float
    baseline: float = 0.0
    amplitude_error: float | None = None
    centroid_error: float | None = None
    sigma_error: float | None = None
    baseline_error: float | None = None
    integrated_intensity_error: float | None = None
    x_unit: str | None = None
    intensity_unit: str | None = None
    integrated_intensity_unit: str | None = None

    @property
    def fwhm(self) -> float:
        """Full width at half maximum for the fitted Gaussian."""
        return 2.0 * np.sqrt(2.0 * np.log(2.0)) * self.sigma

    def as_dict(self) -> dict[str, float | str | None]:
        """Return fit parameters and unit metadata for tables or CSV export."""
        return {
            "amplitude": self.amplitude,
            "centroid": self.centroid,
            "sigma": self.sigma,
            "fwhm": self.fwhm,
            "integrated_intensity": self.integrated_intensity,
            "baseline": self.baseline,
            "amplitude_error": self.amplitude_error,
            "centroid_error": self.centroid_error,
            "sigma_error": self.sigma_error,
            "fwhm_error": (
                2.0 * np.sqrt(2.0 * np.log(2.0)) * self.sigma_error
                if self.sigma_error is not None
                else None
            ),
            "baseline_error": self.baseline_error,
            "integrated_intensity_error": self.integrated_intensity_error,
            "x_unit": self.x_unit,
            "intensity_unit": self.intensity_unit,
            "integrated_intensity_unit": self.integrated_intensity_unit,
        }


def gaussian_model(x, amplitude, centroid, sigma, baseline):
    """Gaussian line profile plus constant baseline."""
    return amplitude * np.exp(-0.5 * ((x - centroid) / sigma) ** 2) + baseline


def integrated_intensity_from_gaussian(amplitude: float, sigma: float) -> float:
    """Return analytic Gaussian integrated intensity."""
    return float(amplitude * abs(sigma) * np.sqrt(2.0 * np.pi))


def _unit_to_string(values) -> str | None:
    """Return a compact unit string for an Astropy Quantity-like input."""
    if isinstance(values, u.Quantity):
        return str(values.unit)
    return None


def _as_value_array(values, name: str) -> np.ndarray:
    """Convert numeric or Quantity input to a one-dimensional float array."""
    if isinstance(values, u.Quantity):
        values = values.to_value(values.unit)
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional.")
    return array


def estimate_initial_gaussian_parameters(x_axis, intensity) -> tuple[float, float, float, float]:
    """Estimate initial ``amplitude, centroid, sigma, baseline`` parameters."""
    x = _as_value_array(x_axis, "x_axis")
    y = _as_value_array(intensity, "intensity")
    _validate_gaussian_inputs(x, y)

    edge_count = max(1, min(5, x.size // 6))
    baseline = float(np.nanmedian(np.concatenate([y[:edge_count], y[-edge_count:]])))
    line_values = y - baseline
    peak_index = int(np.nanargmax(line_values))
    amplitude = float(line_values[peak_index])
    if amplitude <= 0:
        amplitude = float(np.nanmax(y) - np.nanmin(y))
        baseline = float(np.nanmin(y))

    centroid = float(x[peak_index])
    positive = np.clip(line_values, a_min=0.0, a_max=None)
    if np.sum(positive) > 0:
        weighted_centroid = float(np.sum(x * positive) / np.sum(positive))
        variance = float(np.sum(positive * (x - weighted_centroid) ** 2) / np.sum(positive))
        centroid = weighted_centroid
        sigma = np.sqrt(max(variance, 0.0))
    else:
        sigma = 0.0

    x_span = float(np.nanmax(x) - np.nanmin(x))
    channel_width = float(np.nanmedian(np.abs(np.diff(np.sort(x)))))
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = max(x_span / 8.0, channel_width)
    sigma = max(float(sigma), channel_width)
    amplitude = max(float(amplitude), np.finfo(float).eps)
    return amplitude, centroid, sigma, baseline


def _validate_gaussian_inputs(x: np.ndarray, y: np.ndarray) -> None:
    """Validate arrays before Gaussian fitting."""
    if x.shape != y.shape:
        raise ValueError("x_axis and intensity must have the same shape.")
    if x.size < 6:
        raise ValueError("At least six spectral samples are required for fitting.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("x_axis and intensity must contain only finite values.")
    if np.nanmax(x) == np.nanmin(x):
        raise ValueError("x_axis must span a nonzero range.")
    if np.nanmax(y) == np.nanmin(y):
        raise ValueError("intensity must have nonzero variation.")


def fit_single_gaussian(x_axis, intensity) -> GaussianLineFit:
    """Fit one Gaussian emission component plus a constant baseline.

    Parameters
    ----------
    x_axis:
        Spectral coordinate array, usually velocity. Astropy quantities are
        accepted and converted to their current unit before fitting.
    intensity:
        Intensity array. Astropy quantities are accepted and converted to their
        current unit before fitting.

    Returns
    -------
    GaussianLineFit
        Best-fit amplitude, centroid, sigma, baseline, integrated intensity,
        and one-sigma uncertainties when the covariance matrix is available.
    """
    x_unit = _unit_to_string(x_axis)
    intensity_unit = _unit_to_string(intensity)
    integrated_intensity_unit = (
        f"{intensity_unit} {x_unit}"
        if x_unit is not None and intensity_unit is not None
        else None
    )

    x = _as_value_array(x_axis, "x_axis")
    y = _as_value_array(intensity, "intensity")
    _validate_gaussian_inputs(x, y)

    initial = estimate_initial_gaussian_parameters(x, y)
    x_min = float(np.nanmin(x))
    x_max = float(np.nanmax(x))
    min_sigma = float(np.nanmedian(np.abs(np.diff(np.sort(x))))) / 10.0
    max_sigma = max(float(x_max - x_min), min_sigma * 10.0)

    try:
        params, covariance = curve_fit(
            gaussian_model,
            x,
            y,
            p0=initial,
            bounds=(
                [0.0, x_min, min_sigma, -np.inf],
                [np.inf, x_max, max_sigma, np.inf],
            ),
            maxfev=20_000,
        )
    except Exception as exc:
        raise RuntimeError("Single-Gaussian fit did not converge.") from exc

    amplitude, centroid, sigma, baseline = [float(value) for value in params]
    sigma = abs(sigma)
    integrated_intensity = integrated_intensity_from_gaussian(amplitude, sigma)

    errors: list[float | None]
    integrated_error = None
    if covariance is not None and covariance.shape == (4, 4):
        diagonal = np.diag(covariance)
        if np.all(np.isfinite(diagonal)) and np.all(diagonal >= 0):
            errors = [float(value) for value in np.sqrt(diagonal)]
            if amplitude > 0 and sigma > 0:
                relative_variance = (errors[0] / amplitude) ** 2 + (errors[2] / sigma) ** 2
                integrated_error = float(integrated_intensity * np.sqrt(relative_variance))
        else:
            errors = [None, None, None, None]
    else:
        errors = [None, None, None, None]

    return GaussianLineFit(
        amplitude=amplitude,
        centroid=centroid,
        sigma=sigma,
        integrated_intensity=integrated_intensity,
        baseline=baseline,
        amplitude_error=errors[0],
        centroid_error=errors[1],
        sigma_error=errors[2],
        baseline_error=errors[3],
        integrated_intensity_error=integrated_error,
        x_unit=x_unit,
        intensity_unit=intensity_unit,
        integrated_intensity_unit=integrated_intensity_unit,
    )


def fit_multiple_gaussians(x_axis, intensity, n_components: int):
    """Fit multiple Gaussian components to blended or complex spectra."""
    raise NotImplementedError("Multi-Gaussian fitting will be implemented later.")
