"""Spectral extraction helpers for pixels, apertures, and regions."""

from typing import Any

import numpy as np
from astropy import units as u


def extract_pixel_spectrum(cube: Any, x: int, y: int):
    """Extract the spectrum from one spatial pixel.

    Parameters
    ----------
    cube:
        A ``spectral_cube.SpectralCube`` with shape ``(spectral, y, x)``.
    x, y:
        Zero-based pixel coordinates in the spatial plane.

    Returns
    -------
    tuple[astropy.units.Quantity, astropy.units.Quantity]
        Spectral-axis coordinates and intensities for the selected pixel.
    """
    if not isinstance(x, int) or not isinstance(y, int):
        raise TypeError("Pixel coordinates x and y must be integers.")

    _, ny, nx = cube.shape
    if not 0 <= x < nx:
        raise IndexError(f"x pixel {x} is outside valid range [0, {nx - 1}].")
    if not 0 <= y < ny:
        raise IndexError(f"y pixel {y} is outside valid range [0, {ny - 1}].")

    spectral_axis = cube.spectral_axis
    intensity = cube.filled_data[:, y, x]
    return spectral_axis, intensity


def extract_region_spectrum(cube: Any, region: Any, statistic: str = "mean"):
    """Extract a spectrum from a rectangular pixel region.

    Parameters
    ----------
    cube:
        A ``spectral_cube.SpectralCube`` with shape ``(spectral, y, x)``.
    region:
        Dictionary describing an inclusive rectangular pixel region:
        ``{"x_min": int, "x_max": int, "y_min": int, "y_max": int}``.
    statistic:
        Either ``"mean"`` for a region-averaged spectrum or ``"sum"`` for a
        summed spectrum.

    Returns
    -------
    tuple[astropy.units.Quantity, astropy.units.Quantity]
        Spectral-axis coordinates and extracted intensities.
    """
    if statistic not in {"mean", "sum"}:
        raise ValueError("statistic must be either 'mean' or 'sum'.")
    if not isinstance(region, dict):
        raise TypeError("region must be a dictionary of pixel bounds.")

    required_keys = {"x_min", "x_max", "y_min", "y_max"}
    missing_keys = required_keys - set(region)
    if missing_keys:
        missing = ", ".join(sorted(missing_keys))
        raise ValueError(f"region is missing required keys: {missing}")

    bounds = {key: region[key] for key in required_keys}
    if not all(isinstance(value, int) for value in bounds.values()):
        raise TypeError("All rectangular region bounds must be integers.")

    x_min = bounds["x_min"]
    x_max = bounds["x_max"]
    y_min = bounds["y_min"]
    y_max = bounds["y_max"]
    if x_min > x_max or y_min > y_max:
        raise ValueError("Region minimum bounds must be <= maximum bounds.")

    _, ny, nx = cube.shape
    if not 0 <= x_min < nx or not 0 <= x_max < nx:
        raise IndexError(f"x bounds must be within [0, {nx - 1}].")
    if not 0 <= y_min < ny or not 0 <= y_max < ny:
        raise IndexError(f"y bounds must be within [0, {ny - 1}].")

    spectral_axis = cube.spectral_axis
    region_data = cube.filled_data[:, y_min : y_max + 1, x_min : x_max + 1]
    if statistic == "mean":
        intensity = np.nanmean(region_data, axis=(1, 2))
    else:
        intensity = np.nansum(region_data, axis=(1, 2))
    return spectral_axis, intensity


def estimate_rms_noise(
    spectrum: Any,
    line_free_window: tuple[float, float],
    return_quantity: bool = False,
) -> float | u.Quantity:
    """Estimate RMS noise from a line-free spectral interval.

    ``spectrum`` may be the ``(spectral_axis, intensity)`` tuple returned by
    :func:`extract_pixel_spectrum` or :func:`extract_region_spectrum`. Numeric
    ``line_free_window`` endpoints are interpreted in the spectral-axis unit.
    If endpoints are Astropy quantities, they are converted to that unit before
    masking.

    Set ``return_quantity=True`` to return the RMS with the same unit as the
    input intensity. The default remains a plain float for backwards
    compatibility.
    """
    spectral_axis, intensity = spectrum
    if len(line_free_window) != 2:
        raise ValueError("line_free_window must contain exactly two endpoints.")

    lower, upper = line_free_window
    axis_unit = getattr(spectral_axis, "unit", None)
    if axis_unit is not None:
        lower_q = lower if isinstance(lower, u.Quantity) else lower * axis_unit
        upper_q = upper if isinstance(upper, u.Quantity) else upper * axis_unit
        lower_value = lower_q.to_value(axis_unit)
        upper_value = upper_q.to_value(axis_unit)
        axis_values = spectral_axis.to_value(axis_unit)
    else:
        lower_value = float(lower)
        upper_value = float(upper)
        axis_values = np.asarray(spectral_axis)

    if lower_value > upper_value:
        lower_value, upper_value = upper_value, lower_value

    mask = (axis_values >= lower_value) & (axis_values <= upper_value)
    if not np.any(mask):
        raise ValueError("The line-free window does not overlap the spectrum.")

    intensity_unit = getattr(intensity, "unit", None)
    values = intensity[mask]
    if isinstance(values, u.Quantity):
        intensity_unit = values.unit
        values = values.to_value(intensity_unit)
    else:
        values = np.asarray(values)

    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("The line-free window contains no finite intensities.")

    residual = values - np.nanmedian(values)
    rms = float(np.sqrt(np.nanmean(residual**2)))
    if return_quantity and intensity_unit is not None:
        return rms * intensity_unit
    return rms
