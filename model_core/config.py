from collections.abc import Mapping
from numbers import Real

import numpy as np


# Residence time for each compartment (CGC, MGC, TGC, TGN)
TAU = np.array([5.56, 5.56, 5.56, 5.56])

# Donor concentrations (uM)
DONOR_CONC = {
    "UDP-GlcNAc": 9143.0,
    "UDP-Gal": 3810.0,
    "CMP-NeuAc": 2286.0,
    "GDP-Fuc": 5000.0,
    "UDP-GalNAc": 3000.0,
}

# Total enzyme concentrations (uM)
TOTAL_ENZYME_CONC = {
    "ManI": 24.8,
    "ManII": 36.4,
    "FucT": 14.3,
    "GnTI": 29.7,
    "GnTII": 18.2,
    "GnTIII": 111.0,
    "GnTIV": 18.7,
    "GnTV": 21.4,
    "GalT": 87.1,
    "SiaT": 48.4,
}

# Compartment distribution ratio (4 compartments x 10 total enzyme groups)
COMPARTMENT_NAMES = ("CGC", "MGC", "TGC", "TGN")

DIST_MATRIX = np.array([
    [0.05, 0.20, 0.05, 0.15, 0.10, 0.10, 0.10, 0.05, 0.05, 0.05],  # CGC
    [0.15, 0.40, 0.15, 0.30, 0.25, 0.20, 0.25, 0.15, 0.20, 0.15],  # MGC
    [0.40, 0.25, 0.40, 0.35, 0.40, 0.35, 0.30, 0.40, 0.35, 0.40],  # TGC
    [0.40, 0.15, 0.40, 0.20, 0.25, 0.35, 0.35, 0.40, 0.40, 0.40],  # TGN
])

# Total enzyme names corresponding to DIST_MATRIX columns
_BASE_ENZYME_NAMES = [
    "ManI", "ManII", "FucT", "GnTI", "GnTII",
    "GnTIII", "GnTIV", "GnTV", "GalT", "SiaT",
]

# Maps 20 specific enzyme names in the network to 10 total enzyme groups
_ENZYME_GROUP_MAP = {
    "ManI": "ManI", "ManII": "ManII", "FucT": "FucT",
    "GnTI": "GnTI", "GnTII": "GnTII", "GnTIII": "GnTIII",
    "GnTIV": "GnTIV", "GnTV": "GnTV",
    "GnTE_Br1": "GnTIII", "GnTE_Br2": "GnTIII",
    "GnTE_Br3": "GnTIII", "GnTE_Br4": "GnTIII",
    "GalT_Br1": "GalT", "GalT_Br2": "GalT",
    "GalT_Br3": "GalT", "GalT_Br4": "GalT",
    "SiaT_Br1": "SiaT", "SiaT_Br2": "SiaT",
    "SiaT_Br3": "SiaT", "SiaT_Br4": "SiaT",
}


def normalize_dist_matrix(dist_matrix) -> np.ndarray:
    """Return a validated copy whose enzyme columns each sum to one."""
    try:
        matrix = np.asarray(dist_matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Enzyme distribution matrix must contain only numbers") from exc

    expected_shape = (len(COMPARTMENT_NAMES), len(_BASE_ENZYME_NAMES))
    if matrix.shape != expected_shape:
        raise ValueError(
            f"Enzyme distribution matrix must have shape {expected_shape}, got {matrix.shape}"
        )
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Enzyme distribution matrix values must be finite")
    if np.any(matrix < 0):
        raise ValueError("Enzyme distribution matrix values must be non-negative")

    column_totals = matrix.sum(axis=0)
    if np.any(column_totals <= 0):
        raise ValueError("Every enzyme distribution must have a positive total")
    return matrix.copy() / column_totals


def _validate_total_enzyme_conc(total_enzyme_conc) -> dict[str, float]:
    """Validate and copy the ten base-enzyme total concentrations."""
    if not isinstance(total_enzyme_conc, Mapping):
        raise ValueError("Total enzyme concentrations must be a mapping")

    expected = set(_BASE_ENZYME_NAMES)
    supplied = set(total_enzyme_conc)
    missing = sorted(expected - supplied)
    unknown = sorted(supplied - expected)
    if missing:
        raise ValueError(f"Missing total enzyme concentrations: {', '.join(missing)}")
    if unknown:
        raise ValueError(f"Unknown total enzyme concentrations: {', '.join(unknown)}")

    validated = {}
    for name in _BASE_ENZYME_NAMES:
        value = total_enzyme_conc[name]
        if isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError(f"Total enzyme concentration for {name} must be numeric")
        numeric = float(value)
        if not np.isfinite(numeric) or numeric < 0:
            raise ValueError(
                f"Total enzyme concentration for {name} must be finite and non-negative"
            )
        validated[name] = numeric
    return validated


def build_enzyme_dist(total_enzyme_conc=None, dist_matrix=None):
    """
    Computes the concentration of each enzyme in every compartment.
    Returns [{enzyme_name: concentration}, ...], length equals number of compartments.

    Optional parameters are request-local. Defaults preserve the standalone model API,
    while validation guarantees that every enzyme's compartment concentrations sum to
    its supplied total concentration.
    """
    totals = _validate_total_enzyme_conc(
        TOTAL_ENZYME_CONC if total_enzyme_conc is None else total_enzyme_conc
    )
    matrix = normalize_dist_matrix(DIST_MATRIX if dist_matrix is None else dist_matrix)
    n_compartments = len(TAU)
    enzyme_dist = []
    for j in range(n_compartments):
        # First compute compartment concentrations for 10 total enzyme groups
        base_dict = {}
        for idx, base_name in enumerate(_BASE_ENZYME_NAMES):
            base_dict[base_name] = totals[base_name] * matrix[j, idx]

        # Expand to 20 specific enzyme names
        full_dict = {}
        for full_name, base_name in _ENZYME_GROUP_MAP.items():
            full_dict[full_name] = base_dict.get(base_name, 0.0)

        enzyme_dist.append(full_dict)
    return enzyme_dist
