"""Model-facing API helpers for the GlycoGarden Software backend.

This module loads the metabolic model once at import time and exposes
functions to run it on demand for arbitrary promoter strengths and enzyme
compartment distributions.
"""

import sys
from collections.abc import Mapping
from collections import OrderedDict
from copy import deepcopy
from datetime import datetime, timezone
from numbers import Real
from pathlib import Path
from threading import Lock

import numpy as np

_backend_dir = Path(__file__).resolve().parent
_model_dir = (_backend_dir.parent / "model_core").resolve()

# ---------------------------------------------------------------------------
# Load the metabolic model once at import time.
# Building the reaction network is expensive, but solving the steady state is
# the real bottleneck (~2 min per run). Reusing the network saves a few seconds
# per prediction.
# ---------------------------------------------------------------------------
sys.path.insert(0, str(_model_dir))
try:
    from config import (
        COMPARTMENT_NAMES,
        DIST_MATRIX,
        DONOR_CONC,
        TAU,
        TOTAL_ENZYME_CONC,
        build_enzyme_dist,
        normalize_dist_matrix,
    )
    from enzyme import build_enzymes
    from glycoform import Glycoform
    from golgi_model import GolgiModel
    from kinetics import KineticsCalculator
    from reaction import ReactionNetwork
finally:
    sys.path.pop(0)

_BASE_ENZYME_NAMES = list(TOTAL_ENZYME_CONC.keys())

PROMOTER_LEVELS = [
    {"key": "low", "label": "Low", "factor": 0.5},
    {"key": "base", "label": "Base", "factor": 1.0},
    {"key": "elevated", "label": "Elevated", "factor": 1.5},
    {"key": "high", "label": "High", "factor": 2.0},
]

# Capture the original baseline distribution before anything mutates it.
_BASELINE_DIST_MATRIX = normalize_dist_matrix(deepcopy(DIST_MATRIX))
_PREDICTION_LOCK = Lock()
_PREDICTION_CACHE = OrderedDict()
_PREDICTION_CACHE_LIMIT = 32

# ---------------------------------------------------------------------------
# Build the reaction network once.
# ---------------------------------------------------------------------------
_enzymes = build_enzymes()
_network = ReactionNetwork()
_initial_structures = [
    Glycoform(9, 0, 0, 0, 0, 0, 0, 0, 0),
    Glycoform(8, 0, 0, 0, 0, 0, 0, 0, 0),
]
_network.generate_network(_initial_structures, _enzymes)
_n_structs = len(_network.structures)
_kinetics = KineticsCalculator(_network)
_initial_feed = np.zeros(_n_structs)
_initial_feed[:2] = 0.5


def _structure_label(glycoform) -> str:
    fields = [
        ("M", glycoform.man),
        ("F", glycoform.fuc),
        ("G", glycoform.gnb),
        ("b1", glycoform.br1),
        ("b2", glycoform.br2),
        ("b3", glycoform.br3),
        ("b4", glycoform.br4),
        ("Ga", glycoform.gal),
        ("S", glycoform.sia),
    ]
    parts = [f"{prefix}{value}" for prefix, value in fields if value > 0]
    return " ".join(parts) if parts else "Empty"


def _build_structures() -> list[dict]:
    structures = []
    for idx, glycoform in _network.structures.items():
        structures.append({
            "id": int(idx),
            "label": _structure_label(glycoform),
            "man": int(glycoform.man),
            "fuc": int(glycoform.fuc),
            "gnb": int(glycoform.gnb),
            "br1": int(glycoform.br1),
            "br2": int(glycoform.br2),
            "br3": int(glycoform.br3),
            "br4": int(glycoform.br4),
            "gal": int(glycoform.gal),
            "sia": int(glycoform.sia),
        })
    return structures


def _extract_top(final_concs: np.ndarray, top_n: int) -> tuple[list[dict], float]:
    tgn = final_concs[:, -1]
    total = float(tgn.sum())
    if total <= 0:
        return [], total
    normalized = tgn / total
    top_idx = np.argsort(normalized)[-top_n:][::-1]
    top = []
    for rank, idx in enumerate(top_idx, start=1):
        value = float(normalized[idx])
        if value <= 0:
            continue
        top.append({"id": int(idx), "value": round(value, 6), "rank": rank})
    return top, total


def _matrix_to_distribution(matrix: np.ndarray) -> dict[str, dict[str, float]]:
    """Serialize a normalized 4x10 matrix into the public named contract."""
    normalized = normalize_dist_matrix(matrix)
    return {
        enzyme_name: {
            compartment: float(normalized[row_idx, col_idx])
            for row_idx, compartment in enumerate(COMPARTMENT_NAMES)
        }
        for col_idx, enzyme_name in enumerate(_BASE_ENZYME_NAMES)
    }


def _distribution_to_dist_matrix(enzyme_distribution) -> np.ndarray:
    """Strictly validate and normalize a named enzyme distribution payload."""
    if not isinstance(enzyme_distribution, Mapping):
        raise ValueError(
            "enzymeDistribution must be an object containing every configured enzyme"
        )

    expected_enzymes = set(_BASE_ENZYME_NAMES)
    supplied_enzymes = set(enzyme_distribution)
    missing_enzymes = sorted(expected_enzymes - supplied_enzymes)
    unknown_enzymes = sorted(supplied_enzymes - expected_enzymes)
    if missing_enzymes:
        raise ValueError(f"Missing enzyme distributions: {', '.join(missing_enzymes)}")
    if unknown_enzymes:
        raise ValueError(f"Unknown enzyme distributions: {', '.join(unknown_enzymes)}")

    matrix = np.zeros((len(COMPARTMENT_NAMES), len(_BASE_ENZYME_NAMES)), dtype=float)
    expected_compartments = set(COMPARTMENT_NAMES)
    for col_idx, enzyme_name in enumerate(_BASE_ENZYME_NAMES):
        distribution = enzyme_distribution[enzyme_name]
        if not isinstance(distribution, Mapping):
            raise ValueError(f"Distribution for {enzyme_name} must be an object")

        supplied_compartments = set(distribution)
        missing = sorted(expected_compartments - supplied_compartments)
        unknown = sorted(supplied_compartments - expected_compartments)
        if missing:
            raise ValueError(
                f"Distribution for {enzyme_name} is missing compartments: {', '.join(missing)}"
            )
        if unknown:
            raise ValueError(
                f"Distribution for {enzyme_name} has unknown compartments: {', '.join(unknown)}"
            )

        for row_idx, compartment in enumerate(COMPARTMENT_NAMES):
            value = distribution[compartment]
            if isinstance(value, bool) or not isinstance(value, Real):
                raise ValueError(
                    f"Distribution value for {enzyme_name}.{compartment} must be numeric"
                )
            numeric = float(value)
            if not np.isfinite(numeric):
                raise ValueError(
                    f"Distribution value for {enzyme_name}.{compartment} must be finite"
                )
            if numeric < 0:
                raise ValueError(
                    f"Distribution value for {enzyme_name}.{compartment} must be non-negative"
                )
            matrix[row_idx, col_idx] = numeric

    return normalize_dist_matrix(matrix)


def _shift_enzyme_mass(matrix: np.ndarray, from_compartments: list[int], to_compartments: list[int], fraction: float) -> np.ndarray:
    """Shift a fraction of enzyme mass between compartment groups."""
    m = matrix.copy().astype(float)
    n_cols = m.shape[1]
    for col in range(n_cols):
        from_total = sum(m[row, col] for row in from_compartments)
        if from_total == 0:
            continue
        shift = from_total * fraction
        for row in from_compartments:
            m[row, col] -= shift * (m[row, col] / from_total)
        to_total = sum(m[row, col] for row in to_compartments)
        if to_total == 0:
            share = shift / len(to_compartments)
            for row in to_compartments:
                m[row, col] += share
        else:
            for row in to_compartments:
                m[row, col] += shift * (m[row, col] / to_total)
    return normalize_dist_matrix(m)


def get_enzyme_presets() -> list[dict]:
    """Return named presets as complete, lossless four-compartment profiles."""
    baseline = _BASELINE_DIST_MATRIX.copy()
    cgc_biased = _shift_enzyme_mass(baseline, [2, 3], [0], 0.20)
    tgn_biased = _shift_enzyme_mass(baseline, [0, 1], [3], 0.20)
    uniform = np.full_like(baseline, 0.25)
    return [
        {
            "key": "baseline",
            "label": "Baseline",
            "description": "Default compartment distribution from the model config.",
            "distribution": _matrix_to_distribution(baseline),
        },
        {
            "key": "cgc-biased",
            "label": "CGC-biased",
            "description": "Shift 20% of late-compartment enzyme mass toward the cis-Golgi.",
            "distribution": _matrix_to_distribution(cgc_biased),
        },
        {
            "key": "tgn-biased",
            "label": "TGN-biased",
            "description": "Shift 20% of early-compartment enzyme mass toward the trans-Golgi network.",
            "distribution": _matrix_to_distribution(tgn_biased),
        },
        {
            "key": "uniform",
            "label": "Uniform",
            "description": "Each compartment receives an equal 0.25 share of every enzyme group.",
            "distribution": _matrix_to_distribution(uniform),
        },
    ]


def _run_model(promoter_factor: float, dist_matrix: np.ndarray, top_n: int = 15) -> dict:
    """Run one steady-state prediction and return the top glycoforms."""
    scaled_enzyme_conc = {name: value * promoter_factor for name, value in TOTAL_ENZYME_CONC.items()}
    enzyme_dist = build_enzyme_dist(scaled_enzyme_conc, dist_matrix)
    model = GolgiModel(_network, _kinetics, TAU, enzyme_dist, DONOR_CONC, _initial_feed)
    final_concs = model.solve_damped(tol=1e-6, max_iter=50, damping=0.3, verbose=False)

    top, total = _extract_top(final_concs, top_n)
    return {
        "total": round(total, 6),
        "top": top,
    }


def _run_model_queued(promoter_factor: float, dist_matrix: np.ndarray, top_n: int = 15) -> tuple[dict, bool]:
    """Serialize expensive solves and reuse a bounded set of identical results."""
    normalized = normalize_dist_matrix(dist_matrix)
    key = (float(promoter_factor), int(top_n), normalized.astype(np.float64).tobytes())
    with _PREDICTION_LOCK:
        if key in _PREDICTION_CACHE:
            cached = _PREDICTION_CACHE.pop(key)
            _PREDICTION_CACHE[key] = cached
            return deepcopy(cached), True

        result = _run_model(promoter_factor, normalized, top_n)
        _PREDICTION_CACHE[key] = deepcopy(result)
        while len(_PREDICTION_CACHE) > _PREDICTION_CACHE_LIMIT:
            _PREDICTION_CACHE.popitem(last=False)
        return result, False


def predict(promoter_key: str, enzyme_distribution, top_n: int = 15) -> dict:
    """
    Run a real-time prediction for the given promoter level and enzyme distribution.

    Parameters
    ----------
    promoter_key : str
        One of 'low', 'base', 'elevated', 'high'.
    enzyme_distribution : dict[str, dict[str, float]]
        Complete named four-compartment distribution for every base enzyme.
    top_n : int
        Number of top glycoforms to return.

    Returns
    -------
    dict
        {'promoter': ..., 'enzymeDistribution': ..., 'total': ..., 'top': [...]}
    """
    promoter = next((p for p in PROMOTER_LEVELS if p["key"] == promoter_key), None)
    if promoter is None:
        raise ValueError(f"Unknown promoter level: {promoter_key}")

    dist_matrix = _distribution_to_dist_matrix(enzyme_distribution)
    result, cache_hit = _run_model_queued(promoter["factor"], dist_matrix, top_n)

    top_ids = {item["id"] for item in result["top"]}
    top_structures = [s for s in _build_structures() if s["id"] in top_ids]

    return {
        "promoter": promoter_key,
        "enzymeDistribution": _matrix_to_distribution(dist_matrix),
        "total": result["total"],
        "top": result["top"],
        "structures": top_structures,
        "cacheHit": cache_hit,
    }


def generate_matrix(top_n: int = 15) -> dict:
    """Generate the full pre-computed glycoform profile matrix (offline use)."""
    enzyme_presets = get_enzyme_presets()
    results = []
    total_runs = len(PROMOTER_LEVELS) * len(enzyme_presets)
    run_idx = 0
    for promoter in PROMOTER_LEVELS:
        for preset in enzyme_presets:
            run_idx += 1
            print(f"[{run_idx}/{total_runs}] Solving promoter={promoter['key']} preset={preset['key']}...")
            run = predict(promoter["key"], preset["distribution"], top_n)
            run["preset"] = preset["key"]
            run.pop("structures", None)
            results.append(run)
            top_value = run["top"][0]["value"] if run["top"] else 0
            print(f"  -> top abundance: {top_value}, total mass: {run['total']:.2f}")

    return {
        "schemaVersion": "2.0.0",
        "generatedAt": datetime.now(timezone.utc).isoformat(),
        "topN": top_n,
        "promoterLevels": PROMOTER_LEVELS,
        "enzymePresets": [{"key": p["key"], "label": p["label"], "description": p["description"]} for p in enzyme_presets],
        "structures": _build_structures(),
        "results": results,
    }
