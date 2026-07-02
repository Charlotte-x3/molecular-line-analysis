# Project Summary: From Molecular Line Datacubes to Physical Conditions

## Motivation

Millimeter and submillimeter molecular-line observations are a standard way to study cold interstellar gas, dense cores, molecular clouds, and star-forming regions. A spectral-line datacube contains both spatial and velocity information, but turning that cube into physical conditions requires a reproducible analysis chain: inspect the FITS/WCS metadata, extract spectra, identify transitions, measure line intensities, and interpret the results with excitation models.

This project is a compact astronomy research portfolio workflow showing that chain in Python. It uses synthetic demonstration data so the repository can be run without internet access or large telescope files, while keeping the analysis structure close to what would be used on real public molecular-line observations.

## Workflow

The completed workflow currently runs from datacube handling through LTE rotational-diagram analysis:

1. Generate a small synthetic FITS molecular-line datacube with a Gaussian emission line, velocity gradient, beam metadata, and radio-velocity WCS.
2. Load the cube with `spectral-cube`, summarize metadata, extract pixel and rectangular-region spectra, estimate RMS noise, and create moment maps.
3. Load a local molecular transition catalog and match the synthetic cube rest frequency to CO J=1-0.
4. Fit a single Gaussian line profile with a constant baseline and record amplitude, centroid velocity, linewidth, FWHM, integrated intensity, uncertainties, and unit metadata.
5. Build an LTE rotational diagram from a separate same-species synthetic multi-transition table, convert integrated intensities into upper-state column densities, and fit the rotational temperature.

## Key Synthetic Results

The LTE rotational-diagram demonstration uses six transitions of a clearly labelled synthetic molecule, `DemoMol`. The integrated intensities were generated from a known optically thin LTE model with an input rotational temperature near 38 K. The current workflow recovers `T_rot = 37.14 K`, close to the synthetic input value, and reports the fitted intercept as `N_tot / Q(T_rot)` rather than claiming an absolute total column density without a real partition function.

Generated products include:

- `figures/synthetic_moment_maps.png`
- `figures/synthetic_gaussian_line_fit.png`
- `figures/synthetic_rotational_diagram.png`
- `results/synthetic_population_diagram.csv`
- `results/synthetic_rotational_diagram_fit_summary.csv`

## Assumptions And Limitations

The rotational-diagram step assumes LTE excitation, optically thin emission, a single temperature component, beam filling factor approximately one, and multiple transitions from the same species and emitting region. Different molecules such as CO, HCN, and HCO+ are not mixed in one rotational diagram.

The current synthetic workflow does not model optical depth corrections, beam dilution, line blending, calibration systematics, multiple temperature components, or non-LTE excitation. These effects can bias line ratios and make a simple straight-line rotational diagram physically misleading.

## Future Extension

The next planned milestone is a RADEX-style non-LTE grid-comparison interface. That extension will compare observed or synthetic line intensities and ratios against model grids in kinetic temperature, gas density, and column density. This will connect the LTE analysis to a more realistic excitation framework for subthermal gas and density-sensitive molecular tracers.
