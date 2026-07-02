"""Moment map generation for molecular-line datacubes."""

from typing import Any

from astropy import units as u


def _spectral_slab(cube: Any, spectral_window: tuple[float, float]):
    """Return a spectral slab using the cube spectral-axis unit."""
    if len(spectral_window) != 2:
        raise ValueError("spectral_window must contain exactly two endpoints.")

    lower, upper = spectral_window
    axis_unit = cube.spectral_axis.unit
    lower_q = lower if isinstance(lower, u.Quantity) else lower * axis_unit
    upper_q = upper if isinstance(upper, u.Quantity) else upper * axis_unit
    if lower_q > upper_q:
        lower_q, upper_q = upper_q, lower_q
    return cube.spectral_slab(lower_q, upper_q)


def make_moment1(cube: Any, spectral_window: tuple[float, float]) -> Any:
    """Create an intensity-weighted velocity map.

    Moment 1 traces the centroid velocity of the line emission at each spatial
    pixel. In a molecular cloud, coherent gradients in this map can indicate
    rotation, shear, expansion, collapse, or outflow kinematics.
    """
    return _spectral_slab(cube, spectral_window).moment(order=1)


def make_moment0(cube: Any, spectral_window: tuple[float, float]) -> Any:
    """Create an integrated-intensity map over a spectral window.

    Moment 0 integrates line brightness over the selected spectral range. For
    optically thin emission and fixed excitation assumptions, it is often used
    as a proxy for molecular column density or gas mass.
    """
    return _spectral_slab(cube, spectral_window).moment(order=0)


def moment0_to_k_kms(moment0: Any) -> Any:
    """Convert a moment 0 map to the common ``K km/s`` display unit.

    ``spectral-cube`` correctly integrates over the native spectral-axis unit.
    For the synthetic demo cube that means ``K m/s`` internally, while radio
    astronomy figures usually label integrated intensity as ``K km/s``. This
    helper keeps the calculation untouched and only converts compatible output
    objects for plotting or reporting.
    """
    target_unit = u.K * u.km / u.s
    try:
        return moment0.to(target_unit)
    except Exception as exc:
        raise u.UnitConversionError(
            f"Cannot convert moment 0 unit {getattr(moment0, 'unit', None)} "
            f"to {target_unit}."
        ) from exc


def make_moment2(cube: Any, spectral_window: tuple[float, float]) -> Any:
    """Create a velocity-dispersion map over a spectral window.

    This returns the square root of the second central moment, commonly called
    the line-width sigma map. It traces thermal, turbulent, and unresolved
    kinematic broadening in the gas.
    """
    return _spectral_slab(cube, spectral_window).linewidth_sigma()
