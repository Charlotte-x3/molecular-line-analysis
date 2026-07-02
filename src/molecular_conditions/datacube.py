"""FITS datacube loading and metadata inspection."""

from pathlib import Path
from typing import Any

import numpy as np
from spectral_cube import SpectralCube


def load_fits_cube(path: str | Path) -> SpectralCube:
    """Load a molecular-line FITS datacube.

    Parameters
    ----------
    path:
        Path to a local FITS cube, typically stored under ``data/raw``.

    Returns
    -------
    spectral_cube.SpectralCube
        Loaded datacube with WCS and unit metadata.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FITS cube does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Expected a FITS file path, got directory: {path}")

    try:
        return SpectralCube.read(path)
    except Exception as exc:
        raise OSError(f"Could not read FITS cube with spectral-cube: {path}") from exc


def summarize_cube_metadata(cube: Any) -> dict[str, Any]:
    """Summarize shape, units, WCS, beam, and spectral-axis metadata.

    The returned dictionary is intended for quick inspection in notebooks and
    tests, not as a full FITS header replacement.
    """
    try:
        spectral_axis = cube.spectral_axis
    except Exception:
        spectral_axis = None

    has_valid_spectral_axis = False
    spectral_unit = None
    spectral_min = None
    spectral_max = None
    if spectral_axis is not None and len(spectral_axis) > 0:
        finite_values = np.isfinite(spectral_axis.to_value(spectral_axis.unit))
        has_valid_spectral_axis = bool(finite_values.all())
        spectral_unit = str(spectral_axis.unit)
        spectral_min = spectral_axis.min()
        spectral_max = spectral_axis.max()

    wcs = getattr(cube, "wcs", None)
    if wcs is not None:
        wcs_summary = {
            "ctype": list(wcs.wcs.ctype),
            "cunit": [str(unit) for unit in wcs.wcs.cunit],
            "crval": [float(value) for value in wcs.wcs.crval],
            "crpix": [float(value) for value in wcs.wcs.crpix],
            "cdelt": [float(value) for value in wcs.wcs.cdelt],
        }
    else:
        wcs_summary = None

    beam_info = None
    try:
        beam = cube.beam
    except Exception:
        beam = None
    if beam is not None:
        beam_info = {
            "major": beam.major,
            "minor": beam.minor,
            "pa": beam.pa,
        }

    return {
        "shape": tuple(cube.shape),
        "unit": str(cube.unit) if getattr(cube, "unit", None) is not None else None,
        "spectral_axis_unit": spectral_unit,
        "spectral_axis_min": spectral_min,
        "spectral_axis_max": spectral_max,
        "wcs_summary": wcs_summary,
        "beam": beam_info,
        "has_valid_spectral_axis": has_valid_spectral_axis,
    }
