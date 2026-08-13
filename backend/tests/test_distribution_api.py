import json
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

BACKEND_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BACKEND_DIR))

import api
import app as flask_app


class DistributionApiTests(unittest.TestCase):
    def setUp(self):
        self.baseline = api._matrix_to_distribution(api._BASELINE_DIST_MATRIX)
        with api._PREDICTION_LOCK:
            api._PREDICTION_CACHE.clear()

    def test_named_distribution_is_strict_and_normalized(self):
        percentages = {
            enzyme: {compartment: value * 100 for compartment, value in profile.items()}
            for enzyme, profile in self.baseline.items()
        }
        matrix = api._distribution_to_dist_matrix(percentages)
        np.testing.assert_allclose(matrix.sum(axis=0), np.ones(len(api._BASE_ENZYME_NAMES)))
        np.testing.assert_allclose(matrix, api._BASELINE_DIST_MATRIX)

    def test_rejects_malformed_named_distributions(self):
        cases = []

        missing_enzyme = json.loads(json.dumps(self.baseline))
        missing_enzyme.pop("ManI")
        cases.append(missing_enzyme)

        unknown_enzyme = json.loads(json.dumps(self.baseline))
        unknown_enzyme["Unknown"] = unknown_enzyme["ManI"]
        cases.append(unknown_enzyme)

        missing_compartment = json.loads(json.dumps(self.baseline))
        missing_compartment["ManI"].pop("CGC")
        cases.append(missing_compartment)

        unknown_compartment = json.loads(json.dumps(self.baseline))
        unknown_compartment["ManI"]["ER"] = 0
        cases.append(unknown_compartment)

        for value in (-1, float("nan"), float("inf"), True, "0.05"):
            invalid_value = json.loads(json.dumps(self.baseline))
            invalid_value["ManI"]["CGC"] = value
            cases.append(invalid_value)

        all_zero = json.loads(json.dumps(self.baseline))
        all_zero["ManI"] = {name: 0 for name in api.COMPARTMENT_NAMES}
        cases.append(all_zero)

        for distribution in cases:
            with self.subTest(distribution=distribution):
                with self.assertRaises(ValueError):
                    api._distribution_to_dist_matrix(distribution)

    def test_presets_preserve_named_compartment_semantics(self):
        presets = {preset["key"]: preset["distribution"] for preset in api.get_enzyme_presets()}
        for profile in presets["uniform"].values():
            self.assertEqual(profile, {name: 0.25 for name in api.COMPARTMENT_NAMES})

        for enzyme in api._BASE_ENZYME_NAMES:
            baseline = presets["baseline"][enzyme]
            cgc = presets["cgc-biased"][enzyme]
            tgn = presets["tgn-biased"][enzyme]
            self.assertGreater(cgc["CGC"], baseline["CGC"])
            self.assertAlmostEqual(cgc["MGC"], baseline["MGC"])
            self.assertGreater(tgn["TGN"], baseline["TGN"])
            self.assertAlmostEqual(tgn["TGC"], baseline["TGC"])
            self.assertAlmostEqual(sum(cgc.values()), 1)
            self.assertAlmostEqual(sum(tgn.values()), 1)

    def test_predict_echoes_normalized_distribution_without_running_solver(self):
        unnormalized = {
            enzyme: {compartment: value * 10 for compartment, value in profile.items()}
            for enzyme, profile in self.baseline.items()
        }
        with mock.patch.object(api, "_run_model", return_value={"total": 1.0, "top": []}):
            result = api.predict("base", unnormalized)
        np.testing.assert_allclose(
            api._distribution_to_dist_matrix(result["enzymeDistribution"]),
            api._distribution_to_dist_matrix(self.baseline),
        )
        self.assertEqual(result["top"], [])

    def test_identical_predictions_use_bounded_cache(self):
        fake_result = {"total": 1.0, "top": []}
        with mock.patch.object(api, "_run_model", return_value=fake_result) as solver:
            first = api.predict("base", self.baseline)
            second = api.predict("base", self.baseline)
        self.assertFalse(first["cacheHit"])
        self.assertTrue(second["cacheHit"])
        solver.assert_called_once()


class FlaskContractTests(unittest.TestCase):
    def setUp(self):
        flask_app.app.config.update(TESTING=True)
        self.client = flask_app.app.test_client()

    def test_config_exposes_complete_distribution_contract(self):
        response = self.client.get("/api/config")
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(data["compartments"], ["CGC", "MGC", "TGC", "TGN"])
        self.assertEqual(set(data["baselineDistribution"]), set(data["enzymeNames"]))
        self.assertIn("distribution", data["enzymePresets"][0])
        self.assertNotIn("bias", data["enzymePresets"][0])

    def test_frontend_serves_four_compartment_controls(self):
        with self.client.get("/") as response:
            self.assertEqual(response.status_code, 200)
            page = response.get_data(as_text=True)
        self.assertIn("Enzyme compartment distribution", page)
        self.assertIn("config.compartments", page)
        self.assertNotIn("Fine-tune enzyme compartment bias", page)

    def test_legacy_bias_request_is_rejected_before_prediction(self):
        with mock.patch.object(flask_app, "predict") as predictor:
            response = self.client.post(
                "/api/predict",
                json={"promoter": "base", "enzymeBias": {"ManI": 0}},
            )
        self.assertEqual(response.status_code, 400)
        self.assertIn("enzymeDistribution is required", response.get_json()["error"])
        predictor.assert_not_called()

    def test_invalid_distribution_is_rejected_before_solver(self):
        with mock.patch.object(api, "_run_model") as solver:
            response = self.client.post(
                "/api/predict",
                json={"promoter": "base", "enzymeDistribution": {"ManI": {}}},
            )
        self.assertEqual(response.status_code, 400)
        solver.assert_not_called()

    def test_valid_prediction_returns_applied_distribution(self):
        distribution = api._matrix_to_distribution(api._BASELINE_DIST_MATRIX * 100)
        with api._PREDICTION_LOCK:
            api._PREDICTION_CACHE.clear()
        with mock.patch.object(api, "_run_model", return_value={"total": 1.0, "top": []}):
            response = self.client.post(
                "/api/predict",
                json={"promoter": "base", "enzymeDistribution": distribution},
            )
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertIn("enzymeDistribution", data)
        self.assertFalse(data["cacheHit"])
        for profile in data["enzymeDistribution"].values():
            self.assertAlmostEqual(sum(profile.values()), 1)

    def test_missing_offline_matrix_is_reported(self):
        response = self.client.get("/api/matrix")
        self.assertEqual(response.status_code, 404)
        self.assertIn("not found", response.get_json()["error"])


if __name__ == "__main__":
    unittest.main()
