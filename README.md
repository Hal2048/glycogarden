# GlycoDesigner

GlycoDesigner is a interactive simulator, part of Fudan iGEM 2026's GlycoGarden. It aims to explore how Golgi-compartment enzyme localization and cell-physiology parameters (donor concentrations, residence time, cisterna volume, protein production rate) affect predicted N-glycoform profiles.

This repository is self-contained: the browser interface, Flask API, precomputed data, and the model runtime required by the API are all included. The complete research notebooks and broader modeling record live in the separate [Fudan model repository](https://gitlab.igem.org/2026/software/fudan/model).

## Architecture

```text
frontend/       Vue 3 and ECharts browser interface
backend/        Flask API, validation, caching, and static-file server
model_core/     Minimal N-glycosylation model runtime used by the API
data/           Output location for optional precomputed preset results
MODELING.md     Scientific modeling description
```

The frontend calls `/api/config` and `/api/predict`. The backend validates and normalizes every enzyme's four-compartment distribution plus the tunable physiology parameters, then calls `model_core` to solve the steady-state glycosylation model.

## Requirements

- Python 3.10 or newer
- A supported version of NumPy and SciPy
- Node.js is optional and is only needed for the frontend unit test

## Run locally

Clone the repository and run:

```bash
python -m venv .venv
```

Activate the virtual environment:

```bash
# Linux/macOS
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1
```

Install the dependencies and start the application:

```bash
pip install -r backend/requirements.txt
python backend/app.py
```

Open <http://localhost:5000>. No second repository or external model path is required.

## Run with Docker

```bash
docker build -t glycogarden .
docker run --rm -p 5000:5000 glycogarden
```

Then open <http://localhost:5000>.

## Deploy the web application

The repository includes a `render.yaml` Blueprint and a production Gunicorn entrypoint. Connect this public GitLab repository in Render and deploy the Blueprint to publish the complete Flask API, frontend, and model runtime as one web service.

The health check uses `/api/config`. Live simulations can take up to several minutes, so the production server uses one worker and a five-minute request timeout to avoid duplicating the model in memory or terminating valid predictions early.

## Prediction API

`POST /api/predict` requires the complete distribution of every configured enzyme and accepts optional overrides for the physiology parameters. Retrieve the public configuration (defaults and bounds) first with `GET /api/config`.

```json
{
  "enzymeDistribution": {
    "ManI": {"CGC": 0.05, "MGC": 0.15, "TGC": 0.40, "TGN": 0.40}
  },
  "donorConcs": {"UDP-GlcNAc": 9200, "UDP-Gal": 3800, "CMP-NeuAc": 2400, "GDP-Fuc": 5000, "GDP-Man": 2000},
  "tau": 5.56,
  "compartmentVolume": 2.5,
  "proteinProdRate": 1000
}
```

The abbreviated example shows one enzyme. Real requests must include every enzyme returned by `/api/config`. Each profile is normalized to 100 percent by both the API adapter and model layer. All physiology fields are optional; defaults come from `model_core/config.py`. The response echoes the applied parameters together with the derived `totGlycanConc = proteinProdRate · tau / compartmentVolume`.

## Tests

From the repository root:

```bash
python -m unittest discover model_core/tests -v
python -m unittest discover backend/tests -v
node frontend/tests/distribution.test.js
```

## Related repository

The [Fudan model repository](https://gitlab.igem.org/2026/software/fudan/model) contains the complete scientific modeling work, including the code, research outputs, references not required to run GlycoGarden. `model_core/` here is the reviewed runtime subset needed to make this software independently reproducible.

## License

Licensed under the [Apache License 2.0](LICENSE).
