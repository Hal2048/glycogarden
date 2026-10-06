"""Model-facing API helpers for the GlycoGarden Software backend.

This module loads the metabolic model once at import time and exposes
functions to run it on demand for arbitrary enzyme compartment distributions
and cell-physiology parameters (donor concentrations, residence time,
cisterna volume, protein production rate).
"""

import ast
import sys
from collections.abc import Mapping
from collections import OrderedDict
from copy import deepcopy
from numbers import Real
from pathlib import Path
from threading import Lock

import numpy as np

_backend_dir = Path(__file__).resolve().parent
_model_dir = (_backend_dir.parent / "model_core").resolve()

# ---------------------------------------------------------------------------
# Load the metabolic model once at import time.
# Building the reaction network is expensive, but solving the steady state is
# the real bottleneck. Reusing the network saves a few seconds per prediction.
# ---------------------------------------------------------------------------
sys.path.insert(0, str(_model_dir))
try:
    from config import (
        COMPARTMENT_NAMES,
        DIST_MATRIX,
        DONOR_CONC,
        ENZYME_RULES,
        PROTEIN_PROD_RATE,
        COMPARTMENT_VOLUME_SINGLE,
        TAU,
        _BASE_ENZYME_NAMES,
        build_enzyme_dist,
        normalize_dist_matrix,
    )
    from enzyme import build_enzymes
    from glycoform import Glycoform
    from golgi_model import GolgiModel
    from iupac import krambeck_to_iupac
    from kinetics import KineticsCalculator
    from reaction import ReactionNetwork
finally:
    sys.path.pop(0)

_BASE_ENZYME_NAMES = list(_BASE_ENZYME_NAMES)
_N_COMPARTMENTS = len(COMPARTMENT_NAMES)

# ---------------------------------------------------------------------------
# Default physiology parameters (surfaced to the UI as editable controls).
# ---------------------------------------------------------------------------
DEFAULT_TAU = float(TAU[0])
DEFAULT_COMPARTMENT_VOLUME = float(COMPARTMENT_VOLUME_SINGLE)
DEFAULT_PROTEIN_PROD_RATE = float(PROTEIN_PROD_RATE)
DEFAULT_DONOR_CONC = {name: float(value) for name, value in DONOR_CONC.items()}

# Donors the user can tune. "H20" is a unit placeholder for hydrolysis
# reactions (ManI/ManII) and stays fixed at 1.
TUNABLE_DONORS = ("UDP-GlcNAc", "UDP-Gal", "CMP-NeuAc", "GDP-Fuc", "GDP-Man")

# Sensible UI bounds for each tunable parameter.
DONOR_BOUNDS = {
    "UDP-GlcNAc": (100.0, 20000.0),
    "UDP-Gal":    (100.0, 10000.0),
    "CMP-NeuAc":  (50.0,  5000.0),
    "GDP-Fuc":    (100.0, 10000.0),
    "GDP-Man":    (100.0, 10000.0),
}
TAU_BOUNDS = (0.1, 30.0)
COMPARTMENT_VOLUME_BOUNDS = (0.1, 20.0)
PROTEIN_PROD_RATE_BOUNDS = (10.0, 10000.0)

# Capture the original baseline distribution before anything mutates it.
_BASELINE_DIST_MATRIX = normalize_dist_matrix(deepcopy(DIST_MATRIX))
_PREDICTION_LOCK = Lock()
_PREDICTION_CACHE = OrderedDict()
_PREDICTION_CACHE_LIMIT = 32

# ---------------------------------------------------------------------------
# Build the reaction network once.
# ---------------------------------------------------------------------------
_enzymes = build_enzymes(ENZYME_RULES)
_network = ReactionNetwork()
_initial_structures = [
    Glycoform(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    Glycoform(2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
]
_network.generate_network(_initial_structures, _enzymes)
_n_structs = len(_network.structures)
_kinetics = KineticsCalculator(_network)


def _export_network_structures(path: Path) -> None:
    """Persist the authoritative id -> glycoform table for label lookup.

    The layout mirrors ``export_network`` in ``model_workplace/main.py``:
    one line per structure, ``id<TAB>(12-tuple)``. The IUPAC label table
    used by :func:`_build_structures` is built by reading this file back,
    so the file is the single source of truth for id -> structure mapping.
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("GlycoGarden Network Structures (auto-generated at startup)\n")
        f.write(f"Total Structures: {len(_network.structures)}\n")
        f.write("=" * 80 + "\n\n")
        for idx, glycoform in _network.structures.items():
            f.write(f"{idx}\t{tuple(glycoform)}\n")


_NETWORK_STRUCTURES_PATH = _model_dir / "network_structures.txt"
_export_network_structures(_NETWORK_STRUCTURES_PATH)

_STRUCTURE_FIELDS = (
    "man1", "man2", "man3", "fuc", "gnb",
    "br1", "br2", "br3", "br4", "gal", "sia", "br_o",
)


def _read_structure_tuples(path: Path) -> dict[int, tuple]:
    """Parse ``network_structures.txt`` into ``{id: 12-tuple}``.

    This file is the single source of truth for the id -> glycoform mapping.
    Every structure record served by the API is built by reading the tuple
    from this file and running ``krambeck_to_iupac`` on it, so the label
    displayed for a given id always reflects exactly what was persisted.
    """
    tuples: dict[int, tuple] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("=") or line.startswith("Total"):
                continue
            struct_id_str, _, tuple_str = line.partition("\t")
            try:
                struct_id = int(struct_id_str)
                glycoform_tuple = ast.literal_eval(tuple_str)
            except (ValueError, SyntaxError):
                continue
            tuples[struct_id] = glycoform_tuple
    return tuples


def _build_structures() -> list[dict]:
    """Build the structure payload straight from ``network_structures.txt``.

    For each id, the 12-tuple is read from the exported file and converted
    with ``krambeck_to_iupac``; the 12 individual fields are also populated
    from that same tuple so nothing comes from the in-memory network.
    """
    tuples = _read_structure_tuples(_NETWORK_STRUCTURES_PATH)
    structures = []
    for struct_id in sorted(tuples):
        tup = tuples[struct_id]
        record = {"id": struct_id, "label": krambeck_to_iupac(tup)}
        record.update({field: int(value) for field, value in zip(_STRUCTURE_FIELDS, tup)})
        structures.append(record)
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


def _validate_donor_concs(donor_concs) -> dict[str, float]:
    """Validate user-supplied donor concentrations and merge over defaults."""
    if donor_concs is None:
        return {name: DEFAULT_DONOR_CONC[name] for name in TUNABLE_DONORS}
    if not isinstance(donor_concs, Mapping):
        raise ValueError("donorConcs must be an object mapping donor names to concentrations")

    supplied = set(donor_concs)
    unknown = sorted(supplied - set(TUNABLE_DONORS))
    if unknown:
        raise ValueError(f"Unknown donor names in donorConcs: {', '.join(unknown)}")

    merged = {name: DEFAULT_DONOR_CONC[name] for name in TUNABLE_DONORS}
    for name in TUNABLE_DONORS:
        if name not in donor_concs:
            continue
        value = donor_concs[name]
        if isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError(f"Donor concentration for {name} must be numeric")
        numeric = float(value)
        if not np.isfinite(numeric) or numeric <= 0:
            raise ValueError(
                f"Donor concentration for {name} must be finite and positive"
            )
        merged[name] = numeric
    return merged


def _validate_scalar(value, name: str, bounds, default: float) -> float:
    """Validate an optional positive scalar; fall back to the default."""
    if value is None:
        return default
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be numeric")
    numeric = float(value)
    if not np.isfinite(numeric) or numeric <= 0:
        raise ValueError(f"{name} must be finite and positive")
    lo, hi = bounds
    if numeric < lo or numeric > hi:
        raise ValueError(f"{name} must be within [{lo}, {hi}]")
    return numeric


def _run_model(
    dist_matrix: np.ndarray,
    donor_concs: dict[str, float],
    tau_scalar: float,
    compartment_volume: float,
    protein_prod_rate: float,
    top_n: int = 15,
) -> dict:
    """Run one steady-state prediction and return the top glycoforms."""
    tau = np.full(_N_COMPARTMENTS, tau_scalar)
    enzyme_dist = build_enzyme_dist(dist_matrix)

    # Build the full donor table (including the fixed H20 placeholder).
    full_donor_conc = dict(DEFAULT_DONOR_CONC)
    full_donor_conc.update(donor_concs)

    # Total glycan concentration per compartment from physiology parameters
    # (Krambeck-style): influx glycan mass = q_p * tau / V_cisterna.
    tot_glycan_conc = protein_prod_rate * tau_scalar / compartment_volume
    initial_feed = np.zeros(_n_structs)
    initial_feed[:2] = 0.5 * tot_glycan_conc

    model = GolgiModel(_network, _kinetics, tau, enzyme_dist, full_donor_conc, initial_feed)
    final_concs = model.solve_sparse(tol=1e-6, max_iter=50, verbose=False)

    top, total = _extract_top(final_concs, top_n)
    return {
        "total": round(total, 6),
        "top": top,
    }


def _run_model_queued(
    dist_matrix: np.ndarray,
    donor_concs: dict[str, float],
    tau_scalar: float,
    compartment_volume: float,
    protein_prod_rate: float,
    top_n: int = 15,
) -> tuple[dict, bool]:
    """Serialize expensive solves and reuse a bounded set of identical results."""
    normalized = normalize_dist_matrix(dist_matrix)
    donor_key = tuple(sorted(donor_concs.items()))
    key = (
        normalized.astype(np.float64).tobytes(),
        donor_key,
        float(tau_scalar),
        float(compartment_volume),
        float(protein_prod_rate),
        int(top_n),
    )
    with _PREDICTION_LOCK:
        if key in _PREDICTION_CACHE:
            cached = _PREDICTION_CACHE.pop(key)
            _PREDICTION_CACHE[key] = cached
            return deepcopy(cached), True

        result = _run_model(
            normalized, donor_concs, tau_scalar, compartment_volume, protein_prod_rate, top_n
        )
        _PREDICTION_CACHE[key] = deepcopy(result)
        while len(_PREDICTION_CACHE) > _PREDICTION_CACHE_LIMIT:
            _PREDICTION_CACHE.popitem(last=False)
        return result, False


def predict(
    enzyme_distribution,
    donor_concs=None,
    tau=None,
    compartment_volume=None,
    protein_prod_rate=None,
    top_n: int = 15,
) -> dict:
    """
    Run a real-time prediction for the given enzyme distribution and
    cell-physiology parameters.

    Parameters
    ----------
    enzyme_distribution : dict[str, dict[str, float]]
        Complete named four-compartment distribution for every base enzyme.
    donor_concs : dict[str, float], optional
        Override defaults for any of the five tunable donors.
    tau, compartment_volume, protein_prod_rate : float, optional
        Residence time per cisterna (min), single-cisterna volume (uL), and
        protein production rate (uM/min). Defaults come from config.py.
    top_n : int
        Number of top glycoforms to return.

    Returns
    -------
    dict
        The applied parameters plus ``total`` and ``top`` glycoforms.
    """
    dist_matrix = _distribution_to_dist_matrix(enzyme_distribution)
    donor = _validate_donor_concs(donor_concs)
    tau_v = _validate_scalar(tau, "tau", TAU_BOUNDS, DEFAULT_TAU)
    vol_v = _validate_scalar(
        compartment_volume, "compartmentVolume", COMPARTMENT_VOLUME_BOUNDS,
        DEFAULT_COMPARTMENT_VOLUME,
    )
    rate_v = _validate_scalar(
        protein_prod_rate, "proteinProdRate", PROTEIN_PROD_RATE_BOUNDS,
        DEFAULT_PROTEIN_PROD_RATE,
    )

    result, cache_hit = _run_model_queued(dist_matrix, donor, tau_v, vol_v, rate_v, top_n)

    top_ids = {item["id"] for item in result["top"]}
    top_structures = [s for s in _build_structures() if s["id"] in top_ids]

    return {
        "enzymeDistribution": _matrix_to_distribution(dist_matrix),
        "donorConcs": donor,
        "tau": tau_v,
        "compartmentVolume": vol_v,
        "proteinProdRate": rate_v,
        "totGlycanConc": round(rate_v * tau_v / vol_v, 6),
        "total": result["total"],
        "top": result["top"],
        "structures": top_structures,
        "cacheHit": cache_hit,
    }


def get_config_payload() -> dict:
    """Return the dashboard contract: defaults and bounds for every tunable."""
    return {
        "enzymeNames": list(_BASE_ENZYME_NAMES),
        "compartments": list(COMPARTMENT_NAMES),
        "baselineDistribution": _matrix_to_distribution(_BASELINE_DIST_MATRIX),
        "donorNames": list(TUNABLE_DONORS),
        "donorDefaults": {name: DEFAULT_DONOR_CONC[name] for name in TUNABLE_DONORS},
        "donorBounds": {name: list(bounds) for name, bounds in DONOR_BOUNDS.items()},
        "tauDefault": DEFAULT_TAU,
        "tauBounds": list(TAU_BOUNDS),
        "compartmentVolumeDefault": DEFAULT_COMPARTMENT_VOLUME,
        "compartmentVolumeBounds": list(COMPARTMENT_VOLUME_BOUNDS),
        "proteinProdRateDefault": DEFAULT_PROTEIN_PROD_RATE,
        "proteinProdRateBounds": list(PROTEIN_PROD_RATE_BOUNDS),
        "totGlycanConcDefault": round(
            DEFAULT_PROTEIN_PROD_RATE * DEFAULT_TAU / DEFAULT_COMPARTMENT_VOLUME, 6
        ),
    }
