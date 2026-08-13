#!/usr/bin/env python3
"""
Generate the glycoform profile matrix for the GlycoGarden Software dashboard.

Usage:
    cd backend
    python generate_matrix.py

This writes ../data/matrix.json by calling the model wrapper in api.py.
"""

import json
import sys
import traceback
from pathlib import Path

from api import generate_matrix


def main():
    data_dir = Path(__file__).resolve().parent / ".." / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "matrix.json"

    try:
        matrix = generate_matrix(top_n=15)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(matrix, f, indent=2, ensure_ascii=False)
        print(f"Matrix written to: {output_path}")
        print(f"Runs: {len(matrix['results'])}, Structures: {len(matrix['structures'])}")

        assert len(matrix["results"]) == 16
        for run in matrix["results"]:
            assert len(run["top"]) > 0
            top_sum = sum(item["value"] for item in run["top"])
            assert 0.3 <= top_sum <= 1.0, f"Top abundances sum out of range: {top_sum}"
        print("Sanity checks passed.")

    except Exception as e:
        print("\nERROR: failed to generate glycoform matrix.")
        print(str(e))
        traceback.print_exc()
        diagnostic_path = data_dir / "matrix-diagnostic.json"
        with open(diagnostic_path, "w", encoding="utf-8") as f:
            json.dump({
                "schemaVersion": "diagnostic",
                "generatedAt": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
                "error": str(e),
                "details": traceback.format_exc(),
            }, f, indent=2, ensure_ascii=False)
        print(f"\nDiagnostic written to: {diagnostic_path}")
        sys.exit(1)


if __name__ == "__main__":
    main()
