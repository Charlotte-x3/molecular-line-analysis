# Notebooks

The project currently provides three runnable analysis notebooks:

1. `01_datacube_loading_and_moment_maps.ipynb`
   - generates and loads the synthetic FITS cube;
   - inspects metadata, extracts spectra, estimates RMS noise, and creates
     moment maps.
2. `02_line_identification_and_gaussian_fitting.ipynb`
   - matches the synthetic cube rest frequency to a small local catalog;
   - fits a single Gaussian line profile and examines residuals.
3. `03_rotational_diagram_LTE.ipynb`
   - loads a separate synthetic same-species multi-transition table;
   - constructs an LTE population diagram and fits the rotational temperature.

Run the notebooks in this order from the repository root. The Gaussian fit and
LTE rotational diagram currently use two independent synthetic datasets.

A RADEX-style non-LTE notebook is not present because the corresponding model
grid workflow has not been implemented.
