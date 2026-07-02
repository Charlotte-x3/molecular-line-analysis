"""Synthetic FITS datacube generation for runnable examples and tests."""

from pathlib import Path

import numpy as np
from astropy.io import fits


def create_synthetic_line_cube(
    output_path: str | Path,
    overwrite: bool = True,
    seed: int = 42,
) -> Path:
    """Create a small synthetic molecular-line FITS datacube.

    The cube has dimensions ``(velocity, y, x) = (64, 32, 32)``. It contains a
    single Gaussian emission line in brightness-temperature units, a smooth
    spatial amplitude pattern, a linear velocity gradient, and Gaussian random
    noise. The FITS header includes a minimal celestial plus radio-velocity WCS
    that can be read by ``spectral_cube.SpectralCube``.

    Parameters
    ----------
    output_path:
        Destination FITS path.
    overwrite:
        Whether to overwrite an existing file.
    seed:
        Random seed used for reproducible Gaussian noise.

    Returns
    -------
    pathlib.Path
        Path to the written FITS file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    n_velocity, ny, nx = 64, 32, 32

    velocity_m_s = np.linspace(-6_300.0, 6_300.0, n_velocity)
    channel_width = velocity_m_s[1] - velocity_m_s[0]

    y_grid, x_grid = np.mgrid[0:ny, 0:nx]
    x0, y0 = 16.0, 15.0
    spatial_sigma_pix = 6.0
    amplitude_k = 0.35 + 4.8 * np.exp(
        -((x_grid - x0) ** 2 + (y_grid - y0) ** 2) / (2.0 * spatial_sigma_pix**2)
    )

    velocity_gradient_m_s = 2_600.0 * (x_grid - (nx - 1) / 2.0) / (nx - 1)
    line_sigma_m_s = 750.0
    line_profile = amplitude_k[None, :, :] * np.exp(
        -0.5
        * ((velocity_m_s[:, None, None] - velocity_gradient_m_s[None, :, :])
           / line_sigma_m_s)
        ** 2
    )

    noise_rms_k = 0.12
    data = line_profile + rng.normal(0.0, noise_rms_k, size=line_profile.shape)
    data = data.astype("float32")

    header = fits.Header()
    header["SIMPLE"] = True
    header["BITPIX"] = -32
    header["NAXIS"] = 3
    header["NAXIS1"] = nx
    header["NAXIS2"] = ny
    header["NAXIS3"] = n_velocity
    header["BUNIT"] = "K"

    header["CTYPE1"] = "RA---TAN"
    header["CUNIT1"] = "deg"
    header["CRVAL1"] = 83.82208
    header["CRPIX1"] = (nx + 1) / 2.0
    header["CDELT1"] = -1.0 / 3600.0

    header["CTYPE2"] = "DEC--TAN"
    header["CUNIT2"] = "deg"
    header["CRVAL2"] = -5.39111
    header["CRPIX2"] = (ny + 1) / 2.0
    header["CDELT2"] = 1.0 / 3600.0

    header["CTYPE3"] = "VRAD"
    header["CUNIT3"] = "m/s"
    header["CRVAL3"] = float(velocity_m_s[0])
    header["CRPIX3"] = 1.0
    header["CDELT3"] = float(channel_width)
    header["SPECSYS"] = "LSRK"
    header["RESTFRQ"] = 115.2712018e9

    header["BMAJ"] = 3.0 / 3600.0
    header["BMIN"] = 2.2 / 3600.0
    header["BPA"] = 30.0
    header["OBJECT"] = "Synthetic molecular-line cube"
    header["TELESCOP"] = "Synthetic demo"

    fits.PrimaryHDU(data=data, header=header).writeto(output_path, overwrite=overwrite)
    return output_path
