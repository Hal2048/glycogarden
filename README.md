# GlycoGarden

GlycoGarden is Fudan iGEM 2026's interactive simulator for exploring how promoter strength and Golgi-compartment enzyme localization affect predicted N-glycoform profiles.

This repository is self-contained: the browser interface, Flask API, precomputed data, and the model runtime required by the API are all included. The complete research notebooks and broader modeling record live in the separate [Fudan model repository](https://gitlab.igem.org/2026/software/fudan/model).

## Architecture

```text
frontend/       Vue 3 and ECharts browser interface
backend/        Flask API, validation, caching, and static-file server
model_core/     Minimal N-glycosylation model runtime used by the API
data/           Output location for optional precomputed preset results
MODELING.md     Scientific modeling description
```

The frontend calls `/api/config` and `/api/predict`. The backend validates and normalizes every enzyme's four-compartment distribution, then calls `model_core` to solve the steady-state glycosylation model.

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

## Prediction API

`POST /api/predict` requires a promoter level and the complete distribution of every configured enzyme. Retrieve the public configuration first with `GET /api/config`.

```json
{
  "promoter": "base",
  "enzymeDistribution": {
    "ManI": {"CGC": 0.05, "MGC": 0.15, "TGC": 0.40, "TGN": 0.40}
  }
}
```

The abbreviated example shows one enzyme. Real requests must include every enzyme returned by `/api/config`. Each profile is normalized to 100 percent by both the API adapter and model layer.

## Tests

From the repository root:

```bash
python -m unittest discover model_core/tests -v
python -m unittest discover backend/tests -v
node frontend/tests/distribution.test.js
```

## Regenerate preset data

Generating an optional preset/offline matrix for every preset and promoter combination is computationally expensive:

```bash
python backend/generate_matrix.py
```

## Related repository

The [Fudan model repository](https://gitlab.igem.org/2026/software/fudan/model) contains the complete scientific modeling work, including notebooks, research outputs, references, and models not required to run GlycoGarden. `model_core/` here is the reviewed runtime subset needed to make this software independently reproducible.

## License

Licensed under the [Apache License 2.0](LICENSE).
