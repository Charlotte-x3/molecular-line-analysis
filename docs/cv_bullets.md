# CV And Portfolio Text

## Concise CV Bullets

- Built a reproducible Python workflow for molecular-line FITS datacubes, including synthetic data generation, WCS-aware cube loading, spectral extraction, RMS estimation, and moment-map production.
- Implemented molecular transition matching and Gaussian line fitting with uncertainty and unit metadata, producing integrated intensities suitable for downstream excitation analysis.
- Developed an LTE rotational-diagram module for same-species multi-transition measurements, recovering a synthetic input rotational temperature of approximately 38 K under optically thin assumptions.

## Research Portfolio Description

Created a modular Python research portfolio project, "From Molecular Line
Datacubes to Physical Conditions," demonstrating a staged molecular
astrophysics analysis workflow. The project reads a synthetic spectral-line
FITS cube, extracts pixel and region spectra, generates moment maps, identifies
a molecular transition, and fits a Gaussian emission line. A separate
same-species synthetic multi-transition dataset demonstrates LTE
rotational-diagram fitting. The workflow is reproducible without large external
data files and explicitly documents assumptions such as LTE, optically thin
emission, beam filling, and the same-species requirement. RADEX-style non-LTE
modelling remains a planned extension and is not implemented.

## GitHub README Tagline

Reproducible Python workflow from molecular-line FITS datacubes to Gaussian line measurements and LTE rotational-diagram physical conditions.
