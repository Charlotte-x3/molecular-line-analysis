"""Molecular transition catalog handling."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


REQUIRED_TRANSITION_COLUMNS = {
    "species",
    "quantum_numbers",
    "rest_frequency_ghz",
    "upper_energy_k",
    "einstein_a_s",
    "degeneracy_upper",
}


@dataclass(frozen=True)
class MolecularTransition:
    """Minimal metadata for one molecular rotational transition."""

    species: str
    quantum_numbers: str
    rest_frequency_ghz: float
    upper_energy_k: float
    einstein_a_s: float | None = None
    degeneracy_upper: float | None = None


def load_transition_table(path: str | Path) -> pd.DataFrame:
    """Load a local molecular transition catalog.

    Parameters
    ----------
    path:
        CSV file with transition metadata. The expected columns are
        ``species``, ``quantum_numbers``, ``rest_frequency_ghz``,
        ``upper_energy_k``, ``einstein_a_s``, and ``degeneracy_upper``.

    Returns
    -------
    pandas.DataFrame
        Validated transition table sorted by rest frequency.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Transition catalog does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"Expected a transition catalog file, got directory: {path}")

    catalog = pd.read_csv(path)
    missing_columns = REQUIRED_TRANSITION_COLUMNS - set(catalog.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Transition catalog is missing required columns: {missing}")

    catalog = catalog.copy()
    numeric_columns = [
        "rest_frequency_ghz",
        "upper_energy_k",
        "einstein_a_s",
        "degeneracy_upper",
    ]
    for column in numeric_columns:
        catalog[column] = pd.to_numeric(catalog[column], errors="raise")

    catalog["species"] = catalog["species"].astype(str)
    catalog["quantum_numbers"] = catalog["quantum_numbers"].astype(str)
    return catalog.sort_values("rest_frequency_ghz").reset_index(drop=True)


def match_transitions(
    observed_frequency_ghz: float,
    catalog,
    tolerance_mhz: float = 5.0,
) -> pd.DataFrame:
    """Match an observed line frequency to candidate catalog transitions.

    Parameters
    ----------
    observed_frequency_ghz:
        Observed or rest-frame frequency in GHz.
    catalog:
        A transition ``DataFrame`` returned by :func:`load_transition_table`, or
        a list of :class:`MolecularTransition` objects.
    tolerance_mhz:
        Maximum absolute frequency offset in MHz.

    Returns
    -------
    pandas.DataFrame
        Matching rows with additional ``frequency_offset_mhz`` and
        ``abs_frequency_offset_mhz`` columns, sorted by closest match.
    """
    observed_frequency_ghz = float(observed_frequency_ghz)
    tolerance_mhz = float(tolerance_mhz)
    if tolerance_mhz < 0:
        raise ValueError("tolerance_mhz must be non-negative.")

    if isinstance(catalog, pd.DataFrame):
        table = catalog.copy()
    else:
        table = pd.DataFrame(
            [
                {
                    "species": transition.species,
                    "quantum_numbers": transition.quantum_numbers,
                    "rest_frequency_ghz": transition.rest_frequency_ghz,
                    "upper_energy_k": transition.upper_energy_k,
                    "einstein_a_s": transition.einstein_a_s,
                    "degeneracy_upper": transition.degeneracy_upper,
                }
                for transition in catalog
            ]
        )

    missing_columns = REQUIRED_TRANSITION_COLUMNS - set(table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Catalog is missing required columns: {missing}")

    frequency_offset_mhz = (
        observed_frequency_ghz - pd.to_numeric(table["rest_frequency_ghz"])
    ) * 1000.0
    matches = table.copy()
    matches["frequency_offset_mhz"] = frequency_offset_mhz
    matches["abs_frequency_offset_mhz"] = frequency_offset_mhz.abs()
    matches = matches[matches["abs_frequency_offset_mhz"] <= tolerance_mhz]
    return matches.sort_values("abs_frequency_offset_mhz").reset_index(drop=True)
