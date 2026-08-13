import sys
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

MODEL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODEL_DIR))

import config


class DistributionConfigTests(unittest.TestCase):
    def test_normalizes_each_enzyme_column(self):
        matrix = config.DIST_MATRIX * np.arange(1, 11)
        normalized = config.normalize_dist_matrix(matrix)
        np.testing.assert_allclose(normalized.sum(axis=0), np.ones(10))
        np.testing.assert_allclose(normalized, config.DIST_MATRIX)

    def test_build_conserves_supplied_total_concentrations(self):
        totals = {name: value * 1.5 for name, value in config.TOTAL_ENZYME_CONC.items()}
        matrix = config.DIST_MATRIX * 100
        distribution = config.build_enzyme_dist(totals, matrix)
        for name, total in totals.items():
            representative = next(
                full_name
                for full_name, base_name in config._ENZYME_GROUP_MAP.items()
                if base_name == name
            )
            self.assertAlmostEqual(
                sum(compartment[representative] for compartment in distribution),
                total,
            )

    def test_default_call_does_not_mutate_model_configuration(self):
        original_matrix = config.DIST_MATRIX.copy()
        original_totals = config.TOTAL_ENZYME_CONC.copy()
        config.build_enzyme_dist()
        np.testing.assert_array_equal(config.DIST_MATRIX, original_matrix)
        self.assertEqual(config.TOTAL_ENZYME_CONC, original_totals)

    def test_rejects_invalid_matrices(self):
        invalid_matrices = [
            np.ones((3, 10)),
            np.full((4, 10), np.nan),
            np.full((4, 10), np.inf),
            -np.ones((4, 10)),
            np.column_stack([np.zeros(4), np.ones((4, 9))]),
            [["not-a-number"] * 10] * 4,
        ]
        for matrix in invalid_matrices:
            with self.subTest(matrix=np.asarray(matrix).shape):
                with self.assertRaises(ValueError):
                    config.normalize_dist_matrix(matrix)

    def test_rejects_invalid_total_concentrations(self):
        cases = []
        missing = config.TOTAL_ENZYME_CONC.copy()
        missing.pop("ManI")
        cases.append(missing)
        unknown = config.TOTAL_ENZYME_CONC.copy()
        unknown["Unknown"] = 1
        cases.append(unknown)
        for value in (-1, np.nan, np.inf, True, "24.8"):
            invalid = config.TOTAL_ENZYME_CONC.copy()
            invalid["ManI"] = value
            cases.append(invalid)

        for totals in cases:
            with self.subTest(totals=totals):
                with self.assertRaises(ValueError):
                    config.build_enzyme_dist(totals, config.DIST_MATRIX)

    def test_concurrent_builds_are_parameter_isolated(self):
        def build(factor, target_row):
            totals = {name: value * factor for name, value in config.TOTAL_ENZYME_CONC.items()}
            matrix = np.zeros_like(config.DIST_MATRIX)
            matrix[target_row, :] = 1
            compartments = config.build_enzyme_dist(totals, matrix)
            return totals, compartments

        requests = [(0.5, 0), (1.0, 1), (1.5, 2), (2.0, 3)] * 4
        with ThreadPoolExecutor(max_workers=4) as executor:
            results = list(executor.map(lambda args: build(*args), requests))

        for (factor, target_row), (totals, compartments) in zip(requests, results):
            self.assertAlmostEqual(compartments[target_row]["ManI"], totals["ManI"])
            self.assertAlmostEqual(
                sum(compartment["ManI"] for compartment in compartments),
                config.TOTAL_ENZYME_CONC["ManI"] * factor,
            )


if __name__ == "__main__":
    unittest.main()
