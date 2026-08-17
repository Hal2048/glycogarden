from pathlib import Path
import os

from flask import Flask, jsonify, request
from flask_cors import CORS

from api import (
    COMPARTMENT_NAMES,
    PROMOTER_LEVELS,
    _BASELINE_DIST_MATRIX,
    _BASE_ENZYME_NAMES,
    _matrix_to_distribution,
    generate_matrix,
    get_enzyme_presets,
    predict,
)

app = Flask(__name__, static_folder='../frontend', static_url_path='')
CORS(app)

DATA_DIR = Path(__file__).parent.parent / 'data'


@app.route('/api/config', methods=['GET'])
def get_config():
    """Return promoter levels, enzyme names, and presets for the dashboard."""
    return jsonify({
        "promoterLevels": PROMOTER_LEVELS,
        "enzymeNames": _BASE_ENZYME_NAMES,
        "compartments": list(COMPARTMENT_NAMES),
        "baselineDistribution": _matrix_to_distribution(_BASELINE_DIST_MATRIX),
        "enzymePresets": get_enzyme_presets(),
    })


@app.route('/api/predict', methods=['POST'])
def handle_predict():
    """
    Run a real-time prediction.

    Request body:
        {
            "promoter": "base",
            "enzymeDistribution": {
                "ManI": {"CGC": 0.05, "MGC": 0.15, "TGC": 0.4, "TGN": 0.4},
                ...
            }
        }

    Response:
        {
            "promoter": "base",
            "enzymeDistribution": {...},
            "total": 3677492.58,
            "top": [{"id": ..., "value": ..., "rank": 1}, ...]
        }
    """
    data = request.get_json(silent=True) or {}
    promoter = data.get('promoter', 'base')
    if 'enzymeDistribution' not in data:
        return jsonify({
            'error': 'enzymeDistribution is required; the legacy enzymeBias field is no longer supported'
        }), 400
    enzyme_distribution = data.get('enzymeDistribution')

    try:
        result = predict(promoter, enzyme_distribution, top_n=15)
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/api/matrix', methods=['GET'])
def get_matrix():
    """Return the pre-computed glycoform profile matrix (static fallback)."""
    matrix_path = DATA_DIR / 'matrix.json'
    if not matrix_path.exists():
        return jsonify({'error': 'Pre-computed matrix not found'}), 404
    with open(matrix_path, 'r', encoding='utf-8') as f:
        matrix = __import__('json').load(f)
    if matrix.get('schemaVersion') != '2.0.0':
        return jsonify({
            'error': 'Pre-computed matrix schema is outdated; run backend/generate_matrix.py to regenerate schema 2.0.0'
        }), 409
    return jsonify(matrix)


@app.route('/')
def index():
    return app.send_static_file('index.html')


if __name__ == '__main__':
    debug = os.environ.get('GLYCOGARDEN_DEBUG', '').lower() in {'1', 'true', 'yes'}
    port = int(os.environ.get('PORT', '5000'))
    app.run(host='0.0.0.0', port=port, debug=debug)
