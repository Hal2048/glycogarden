from pathlib import Path
import os

from flask import Flask, jsonify, request
from flask_cors import CORS

from api import get_config_payload, predict

app = Flask(__name__, static_folder='../frontend', static_url_path='')
CORS(app)


@app.route('/api/config', methods=['GET'])
def get_config():
    """Return enzyme names, baseline distribution, and physiology defaults."""
    return jsonify(get_config_payload())


@app.route('/api/predict', methods=['POST'])
def handle_predict():
    """
    Run a real-time prediction.

    Request body:
        {
            "enzymeDistribution": {
                "ManI": {"CGC": 0.05, "MGC": 0.15, "TGC": 0.4, "TGN": 0.4},
                ...
            },
            "donorConcs": {"UDP-GlcNAc": 9200, ...},   # optional
            "tau": 5.56,                                # optional
            "compartmentVolume": 2.5,                   # optional
            "proteinProdRate": 1000                     # optional
        }

    Response echoes the applied parameters (including the derived
    ``totGlycanConc``) along with ``total`` and ``top`` glycoforms.
    """
    data = request.get_json(silent=True) or {}
    if 'enzymeDistribution' not in data:
        return jsonify({
            'error': 'enzymeDistribution is required'
        }), 400

    try:
        result = predict(
            data['enzymeDistribution'],
            donor_concs=data.get('donorConcs'),
            enzyme_concs=data.get('enzymeConcs'),
            tau=data.get('tau'),
            compartment_volume=data.get('compartmentVolume'),
            protein_prod_rate=data.get('proteinProdRate'),
            top_n=15,
        )
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/')
def index():
    return app.send_static_file('index.html')


if __name__ == '__main__':
    debug = os.environ.get('GLYCODESIGNER_DEBUG', '').lower() in {'1', 'true', 'yes'}
    port = int(os.environ.get('PORT', '5000'))
    app.run(host='0.0.0.0', port=port, debug=debug)
