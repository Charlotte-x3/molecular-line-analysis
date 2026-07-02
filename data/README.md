# Data Directory

This directory contains the small synthetic datasets used by the runnable
notebooks and reserves separate locations for future local astronomy data.

- `demo/`: committed synthetic demonstration data:
  - `synthetic_line_cube.fits`: a generated CO J=1-0-style FITS cube with WCS,
    beam metadata, a velocity gradient, and Gaussian noise;
  - `transitions_demo.csv`: a small local transition table used for frequency
    matching;
  - `demo_rotational_lines.csv`: six synthetic same-species `DemoMol`
    measurements used by the LTE rotational-diagram notebook.
- `raw/`: local downloaded FITS datacubes or spectra.
- `processed/`: local cleaned spectra, masks, compact tables, or derived
  intermediate products.
- `external/`: local external catalogs or model grids.

The Gaussian-fitting notebook uses the synthetic CO cube, while the LTE
rotational-diagram notebook uses the independent `DemoMol` table. These are two
separate synthetic datasets, not a single connected measurement chain.

Large local astronomy files under `raw/`, `processed/`, and `external/` are
ignored by Git. The repository does not currently download public telescope
data automatically.
