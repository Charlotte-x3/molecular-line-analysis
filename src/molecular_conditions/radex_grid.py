"""RADEX-style non-LTE model-grid comparison interfaces."""

from pathlib import Path


def load_model_grid(path: str | Path):
    """Load a precomputed non-LTE molecular excitation model grid."""
    raise NotImplementedError("Model-grid loading will be implemented later.")


def compare_observed_to_grid(observed_lines, model_grid, metric: str = "chi_square"):
    """Compare observed line intensities or ratios with a model grid."""
    raise NotImplementedError("Grid comparison will be implemented later.")


def summarize_best_fit(comparison_result):
    """Return best-fit kinetic temperature, density, and column density."""
    raise NotImplementedError("Best-fit summary will be implemented later.")
