from collections.abc import Mapping
from numbers import Real

import numpy as np


# ----------------------------------------------------------------------------
# Cell physiology parameters (scalar, user-tunable from the dashboard)
# ----------------------------------------------------------------------------

# Protein production rate (uM / min)
PROTEIN_PROD_RATE = 1000.0

# Volume of a single cisterna (uL)
COMPARTMENT_VOLUME_SINGLE = 2.5

# Residence time per compartment (CGC, MGC, TGC, TGN), in minutes.
tau = 5.56
TAU = np.array([tau, tau, tau, tau])

# Total glycan concentration per compartment (uM), derived from the physiology
# parameters above. The dashboard surfaces the scalar formula
#     TOT_GLYCAN_CONC = PROTEIN_PROD_RATE * tau / COMPARTMENT_VOLUME_SINGLE
# and recomputes it whenever the user edits one of the factors.
TOT_GLYCAN_CONC = PROTEIN_PROD_RATE * TAU / COMPARTMENT_VOLUME_SINGLE  # uM


# ----------------------------------------------------------------------------
# Donor concentrations (uM)
# ----------------------------------------------------------------------------

DONOR_CONC = {
    "UDP-GlcNAc": 9200.0,
    "UDP-Gal": 3800.0,
    "CMP-NeuAc": 2400.0,
    "GDP-Fuc": 5000.0,
    "GDP-Man": 2000.0,
    "H20": 1.0,
}


# ----------------------------------------------------------------------------
# Enzyme kinetic rules and compartment concentrations
# ----------------------------------------------------------------------------

# Total enzyme concentrations (uM) and default compartment distribution ratio
# (CGC, MGC, TGC, TGN), keyed by specific enzyme name.
ENZYME_CONC = {
    "ManI_9":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManI_8":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManI_7":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManI_6":   [1.78, [0.05, 0.15, 0.40, 0.40]],
    "ManII_5":  [1.32, [0.20, 0.40, 0.25, 0.15]],
    "ManII_4":  [1.32, [0.20, 0.40, 0.25, 0.15]],
    "FucT":     [2.50, [0.05, 0.15, 0.40, 0.40]],
    "GnTI":     [3.05, [0.15, 0.30, 0.35, 0.20]],
    "GnTII":    [1.29, [0.10, 0.25, 0.40, 0.25]],
    "GnTIV":    [3.62, [0.10, 0.20, 0.35, 0.35]],
    "GnTV":     [0.4,  [0.10, 0.25, 0.30, 0.35]],
    "GnTE_Br1": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GnTE_Br2": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GnTE_Br3": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GnTE_Br4": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "GalT_Br1": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "GalT_Br2": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "GalT_Br3": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "GalT_Br4": [0.66, [0.05, 0.20, 0.35, 0.40]],
    "SiaT_Br1": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "SiaT_Br2": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "SiaT_Br3": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "SiaT_Br4": [3.47, [0.05, 0.15, 0.40, 0.40]],
    "Och1":     [0,    [0.25, 0.25, 0.25, 0.25]],
    "Mnn9":     [0,    [0.25, 0.25, 0.25, 0.25]],
    "GnTIII":   [0,    [0.25, 0.25, 0.25, 0.25]],
}


ENZYME_RULES = [
    # ManI, 1
    {
        "name": "ManI_M9",
        "condition": "g.man1 == 2 and g.man2 == 2 and g.man3 == 2", # M9
        "product_struc": "Glycoform(g.man1, g.man2, g.man3-1, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManI, 2
    {
        "name": "ManI_M8",
        "condition": "g.man1 == 2 and g.man2 == 2 and g.man3 == 1", # M8
        "product_struc": "Glycoform(g.man1-1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManI, 3
    {
        "name": "ManI_M7",
        "condition": "g.man1 == 1 and g.man2 == 2 and g.man3 == 1", # M7
        "product_struc": "Glycoform(g.man1, g.man2, g.man3-1, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManI, 3
    {
        "name": "ManI_M6",
        "condition": "g.man1 == 1 and g.man2 == 2 and g.man3 == 0", # M6
        "product_struc": "Glycoform(g.man1, g.man2-1, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 100, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManII, 1
    {
        "name": "ManII_M5",
        "condition": "g.man1 ==1 and g.man2 == 1 and g.man3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2-1, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 200, "Kmd": 0, "kf": 1924,
        "adjustments": []
    },

    # ManII, 2
    {
        "name": "ManII_M4",
        "condition": "g.man1 ==1 and g.man2 == 0 and g.man3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1-1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "H20",
        "Km": 200, "Kmd": 0, "kf": 1924,
        "adjustments": [("g.man1 ==1 and g.man2 == 0", 0.50)]
    },

    # FucT
    {
        "name": "FucT",
        "condition": "g.fuc == 0 and g.br4 > 0 and g.gnb == 0 and g.gal == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc+1, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "GDP-Fuc",
        "Km": 25, "Kmd": 46, "kf": 253,
        "adjustments": []
    },

    # GnTI
    {
        "name": "GnTI",
        "condition": "g.br4 == 0 and g.man1 + g.man2 + g.man3 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 260, "Kmd": 170, "kf": 990,
        "adjustments": []
    },

    # GnTII
    {
        "name": "GnTII",
        "condition": "g.br2 == 0 and g.man1 + g.man2 + g.man3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 190, "Kmd": 3100, "kf": 607,
        "adjustments": []
    },

    # # GnTIII
    # {
    #     "name": "GnTIII",
    #     "condition": "g.gnb == 0 and g.br4 > 0 and g.gal == 0",
    #     "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb+1, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
    #     "cosubstrate": "UDP-GlcNAc",
    #     "Km": 4000, "Kmd": 3100, "kf": 607,
    #     "adjustments": [("g.br2 > 0", 0.048)]
    # },

    # GnTIV
    {
        "name": "GnTIV",
        "condition": "g.br3 == 0 and g.br4 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 3400, "Kmd": 8300, "kf": 187,
        "adjustments": [("g.br2 == 0", 5), ("g.br2 > 1 or g.br1 > 1", 1.5), ("g.br1 > 0", 0.178)]
    },

    # GnTV
    {
        "name": "GnTV",
        "condition": "g.br1 == 0 and g.br2 == 1 and g.gnb == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 130, "Kmd": 3500, "kf": 1410,
        "adjustments": [("g.br3 > 0", 0.692)]
    },

    # GnTE, 1
    {
        "name": "GnTE_Br1",
        "condition": "g.br1 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+2, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 700, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br2 == 0", 4)]
    },

    # GnTE, 2
    {
        "name": "GnTE_Br2",
        "condition": "g.br2 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+2, g.br3, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 700, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br1 == 0", 4)]
    },

    # GnTE, 3
    {
        "name": "GnTE_Br3",
        "condition": "g.br3 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+2, g.br4, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 7000, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br4 == 0", 4)]
    },

    # GnTE, 4
    {
        "name": "GnTE_Br4",
        "condition": "g.br4 == 2",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+2, g.gal, g.sia, g.br_o)",
        "cosubstrate": "UDP-GlcNAc",
        "Km": 7000, "Kmd": 55, "kf": 25,
        "adjustments": [("g.br3 == 0", 4)]
    },

    # GalT, 1
    {
        "name": "GalT_Br1",
        "condition": "g.br1 == 1 or g.br1 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 150, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # GalT, 2
    {
        "name": "GalT_Br2",
        "condition": "g.br2 == 1 or g.br2 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 135, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # GalT, 3
    {
        "name": "GalT_Br3",
        "condition": "g.br3 == 1 or g.br3 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 80, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # GalT, 4
    {
        "name": "GalT_Br4",
        "condition": "g.br4 == 1 or g.br4 == 4",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal+1, g.sia, g.br_o)",
        "cosubstrate": "UDP-Gal",
        "Km": 4000, "Kmd": 0, "kf": 8712,
        "adjustments": [("g.gnb > 0 and g.br2 > 0", 3.62)]
    },

    # SiaT, 1
    {
        "name": "SiaT_Br1",
        "condition": "g.br1 == 2 or g.br1 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1+1, g.br2, g.br3, g.br4, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # SiaT, 2
    {
        "name": "SiaT_Br2",
        "condition": "g.br2 == 2 or g.br2 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2+1, g.br3, g.br4, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # SiaT, 3
    {
        "name": "SiaT_Br3",
        "condition": "g.br3 == 2 or g.br3 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3+1, g.br4, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # SiaT, 4
    {
        "name": "SiaT_Br4",
        "condition": "g.br4 == 2 or g.br4 == 5",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4+1, g.gal, g.sia+1, g.br_o)",
        "cosubstrate": "CMP-NeuAc",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },

    # Och1
    {
        "name": "Och1",
        "condition": "g.man1 > 0 and g.man2 > 0 and g.man3 > 0 and g.fuc == 0 and g.gnb == 0 and g.br1 == 0 and g.br2 == 0 and g.br3 == 0 and g.br4 == 0 and g.gal == 0 and g.sia == 0 and g.br_o == 0",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o+1)",
        "cosubstrate": "GDP-Man",
        "Km": 390, "Kmd": 10, "kf": 5.2,
        "adjustments": [("g.man1 + g.man2 + g.man3 == 6", 0.74), ("g.man1 == 2 and g.man2 == 2 and g.man3 == 1", 0.55), ("g.man1 == 1 and g.man2 == 2 and g.man3 == 2", 0.27), ("g.man1 == 2 and g.man2 == 1 and g.man3 == 1", 0.6), ("g.man1 == 1 and g.man2 == 1 and g.man3 == 2", 0.25), ("g.man1 == 1 and g.man2 == 2 and g.man3 == 1", 0.14), ("g.man1 == 1 and g.man2 == 1 and g.man3 == 1", 0.16)]
    },

    # Mnn9
    {
        "name": "Mnn9",
        "condition": "g.br_o == 1",
        "product_struc": "Glycoform(g.man1, g.man2, g.man3, g.fuc, g.gnb, g.br1, g.br2, g.br3, g.br4, g.gal, g.sia, g.br_o+1)",
        "cosubstrate": "GDP-Man",
        "Km": 260, "Kmd": 57, "kf": 484,
        "adjustments": [("g.sia > 1", 5.0)]
    },
]


# ----------------------------------------------------------------------------
# Dashboard compatibility layer
#
# The interactive dashboard exposes twelve base enzyme groups and edits a
# single 4-compartment distribution plus a total concentration per group.
# Each specific enzyme in ENZYME_CONC belongs to exactly one base group and
# shares its distribution profile with the other members of that group, so a
# 4x12 base-enzyme matrix cleanly maps back onto the 23 specific ENZYME_CONC
# entries.
# ----------------------------------------------------------------------------

n_compartments = 4

COMPARTMENT_NAMES = ("CGC", "MGC", "TGC", "TGN")

# Base enzyme groups surfaced in the dashboard distribution editor.
_BASE_ENZYME_NAMES = [
    "ManI", "ManII", "FucT", "GnTI", "GnTII",
    "GnTIII", "GnTIV", "GnTV", "GalT", "SiaT",
    "Och1", "Mnn9",
]

# Maps each specific enzyme name in ENZYME_CONC to a base enzyme group.
_ENZYME_GROUP_MAP = {
    "ManI_9": "ManI", "ManI_8": "ManI", "ManI_7": "ManI", "ManI_6": "ManI",
    "ManII_5": "ManII", "ManII_4": "ManII",
    "FucT": "FucT",
    "GnTI": "GnTI", "GnTII": "GnTII",
    "GnTIV": "GnTIV", "GnTV": "GnTV",
    "GnTE_Br1": "GnTIII", "GnTE_Br2": "GnTIII",
    "GnTE_Br3": "GnTIII", "GnTE_Br4": "GnTIII",
    "GalT_Br1": "GalT", "GalT_Br2": "GalT",
    "GalT_Br3": "GalT", "GalT_Br4": "GalT",
    "SiaT_Br1": "SiaT", "SiaT_Br2": "SiaT",
    "SiaT_Br3": "SiaT", "SiaT_Br4": "SiaT",
    "Och1": "Och1", "Mnn9": "Mnn9",
}

# Specific enzymes with no active rule keep their default distribution.
_UNGROUPPED_ENZYMES = ("GnTIII",)


def _compute_total_enzyme_conc() -> dict:
    """Sum specific-enzyme totals into the ten base groups."""
    totals = {name: 0.0 for name in _BASE_ENZYME_NAMES}
    for specific, (total, _ratios) in ENZYME_CONC.items():
        base = _ENZYME_GROUP_MAP.get(specific)
        if base is not None:
            totals[base] += total
    return totals


TOTAL_ENZYME_CONC = _compute_total_enzyme_conc()


def _compute_baseline_dist_matrix() -> np.ndarray:
    """Derive the baseline 4x10 base-enzyme distribution from ENZYME_CONC."""
    matrix = np.zeros((n_compartments, len(_BASE_ENZYME_NAMES)))
    for k, base in enumerate(_BASE_ENZYME_NAMES):
        representative = next(
            (specific for specific, b in _ENZYME_GROUP_MAP.items() if b == base),
            None,
        )
        if representative is not None:
            matrix[:, k] = ENZYME_CONC[representative][1]
        else:
            matrix[:, k] = 1.0 / n_compartments
    return matrix


DIST_MATRIX = _compute_baseline_dist_matrix()


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


def build_enzyme_dist(dist_matrix=None, concentrations=None):
    """
    Computes the concentration of each specific enzyme in every compartment.

    Parameters
    ----------
    dist_matrix : optional, shape (4, 12)
        Base-enzyme distribution matrix; column ``k`` is the four-compartment
        share of ``_BASE_ENZYME_NAMES[k]``. Each specific enzyme inherits its
        base group's profile. When omitted, the baseline matrix derived from
        ``ENZYME_CONC`` is used.
    concentrations : optional, dict[str, float]
        Total concentration override (uM) per base enzyme group. Members of a
        group are scaled proportionally to their original totals; a group
        whose original total is zero (e.g. Och1, Mnn9) assigns the override
        directly to each of its members.

    Returns
    -------
    list[dict]
        ``[{specific_enzyme_name: concentration}, ...]`` with one dict per
        compartment, matching the layout expected by :class:`GolgiModel`.
    """
    if dist_matrix is None:
        base_matrix = DIST_MATRIX
    else:
        base_matrix = normalize_dist_matrix(dist_matrix)

    base_index = {name: k for k, name in enumerate(_BASE_ENZYME_NAMES)}
    base_totals = {name: 0.0 for name in _BASE_ENZYME_NAMES}
    base_members = {name: [] for name in _BASE_ENZYME_NAMES}
    for specific, (total, _ratios) in ENZYME_CONC.items():
        base = _ENZYME_GROUP_MAP.get(specific)
        if base is not None:
            base_totals[base] += total
            base_members[base].append(specific)

    overrides = concentrations or {}
    enzyme_dist = []
    for j in range(n_compartments):
        comp_dict = {}
        for specific, (total, default_ratios) in ENZYME_CONC.items():
            base = _ENZYME_GROUP_MAP.get(specific)
            if base is not None:
                share = base_matrix[j, base_index[base]]
            else:
                share = default_ratios[j]
            specific_total = total
            if base in overrides:
                base_total = base_totals[base]
                if base_total > 0:
                    specific_total = total * (overrides[base] / base_total)
                else:
                    n_members = max(1, len(base_members[base]))
                    specific_total = overrides[base] / n_members
            comp_dict[specific] = specific_total * share
        enzyme_dist.append(comp_dict)
    return enzyme_dist


if __name__ == "__main__":
    enzyme_dist = build_enzyme_dist()
    for i, comp in enumerate(enzyme_dist):
        print(f"Compartment {i+1}:")
        for enzyme, conc in comp.items():
            print(f"  {enzyme}: {conc:.2f}")
