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

    def test_default_call_does_not_mutate_model_configuration(self):
        original_matrix = config.DIST_MATRIX.copy()
        config.build_enzyme_dist()
        np.testing.assert_array_equal(config.DIST_MATRIX, original_matrix)

    def test_baseline_conserves_enzyme_conc_totals(self):
        distribution = config.build_enzyme_dist()
        for base_name in config._BASE_ENZYME_NAMES:
            representative = next(
                specific
                for specific, group in config._ENZYME_GROUP_MAP.items()
                if group == base_name
            )
            expected_total = config.ENZYME_CONC[representative][0]
            observed_total = sum(comp[representative] for comp in distribution)
            self.assertAlmostEqual(observed_total, expected_total)

    def test_custom_matrix_overrides_baseline(self):
        matrix = np.zeros_like(config.DIST_MATRIX)
        matrix[2, :] = 1.0  # everything in TGC
        distribution = config.build_enzyme_dist(matrix)
        for base_name in config._BASE_ENZYME_NAMES:
            representative = next(
                specific
                for specific, group in config._ENZYME_GROUP_MAP.items()
                if group == base_name
            )
            self.assertAlmostEqual(distribution[2][representative],
                                   config.ENZYME_CONC[representative][0])
            for j in (0, 1, 3):
                self.assertAlmostEqual(distribution[j][representative], 0.0)

    def test_ungrouped_enzymes_keep_default_ratios(self):
        distribution = config.build_enzyme_dist()
        for specific in config._UNGROUPPED_ENZYMES:
            for j in range(config.n_compartments):
                self.assertAlmostEqual(
                    distribution[j][specific],
                    config.ENZYME_CONC[specific][0]
                    * config.ENZYME_CONC[specific][1][j],
                )

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

    def test_concurrent_builds_are_parameter_isolated(self):
        def build(target_row):
            matrix = np.zeros_like(config.DIST_MATRIX)
            matrix[target_row, :] = 1
            return config.build_enzyme_dist(matrix)

        with ThreadPoolExecutor(max_workers=4) as executor:
            results = list(executor.map(build, [0, 1, 2, 3] * 4))

        for target_row, compartments in enumerate(results[:4]):
            self.assertAlmostEqual(
                compartments[target_row]["ManI_9"],
                config.ENZYME_CONC["ManI_9"][0],
            )


if __name__ == "__main__":
    unittest.main()
