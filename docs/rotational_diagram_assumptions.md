# LTE Rotational Diagram Assumptions

This project uses a rotational diagram as a transparent LTE demonstration for a synthetic same-species dataset. The current implementation is designed for teaching, portfolio review, and method scaffolding, not for replacing a full molecular excitation analysis.

## Same-Species Requirement

A rotational diagram compares the relative populations of different upper energy levels for one molecule. The fitted line assumes that all plotted points share one total column density scale and one partition function. Mixing different molecules, such as CO, HCN, and HCO+, would combine unrelated abundances, dipole moments, excitation conditions, and partition functions. That would make the inferred temperature physically meaningless.

The demo therefore uses only synthetic `DemoMol` transitions. Real applications should use multiple transitions from the same species and, as much as possible, the same emitting region and beam treatment.

## Why The Slope Gives T_rot

Under optically thin LTE assumptions,

```text
ln(N_u / g_u) = ln(N_tot / Q(T_rot)) - E_u / T_rot
```

where `E_u` is expressed in Kelvin. This is a straight line with slope `-1 / T_rot`. A steeper negative slope corresponds to a lower rotational temperature, while a shallower slope corresponds to a warmer excitation distribution.

## Why The Code Reports N_tot / Q(T_rot)

The intercept is `ln(N_tot / Q(T_rot))`. A true total column density requires a partition function `Q(T_rot)` for the specific molecule and state conventions used in the line catalog. The synthetic `DemoMol` example does not have a real partition function, so the code reports `column_density_over_partition_function_cm2` instead of pretending to know `N_tot`.

For real molecules, the next step would be to obtain a consistent partition function from a trusted molecular database or spectroscopy reference, then compute `N_tot = exp(intercept) * Q(T_rot)`.

## When The Assumptions Fail

Optical depth: If lines are optically thick, the measured intensity no longer scales linearly with upper-state column density. Points can fall below the optically thin expectation, biasing the fit.

Beam dilution: Different transitions may be observed at different frequencies or resolutions. If the emitting region does not fill the beam equally, line ratios can reflect beam effects rather than level populations.

Line blending: Overlapping spectral features can inflate or distort integrated intensities. Blended lines need careful deblending or exclusion from the diagram.

Calibration uncertainty: Absolute flux calibration can dominate the uncertainty budget, especially when combining lines from different instruments, bands, or observing programs.

Non-LTE excitation: If densities are below critical densities or radiative pumping is important, level populations may not follow a single Boltzmann temperature. In that case, a RADEX-style non-LTE comparison is more appropriate.
