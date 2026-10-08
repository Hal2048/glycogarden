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

    def test_predict_echoes_normalized_distribution_without_running_solver(self):
        unnormalized = {
            enzyme: {compartment: value * 10 for compartment, value in profile.items()}
            for enzyme, profile in self.baseline.items()
        }
        with mock.patch.object(api, "_run_model", return_value={"total": 1.0, "top": []}):
            result = api.predict(unnormalized)
        np.testing.assert_allclose(
            api._distribution_to_dist_matrix(result["enzymeDistribution"]),
            api._distribution_to_dist_matrix(self.baseline),
        )
        self.assertEqual(result["top"], [])
        self.assertAlmostEqual(result["tau"], api.DEFAULT_TAU)
        self.assertAlmostEqual(result["compartmentVolume"], api.DEFAULT_COMPARTMENT_VOLUME)
        self.assertAlmostEqual(result["proteinProdRate"], api.DEFAULT_PROTEIN_PROD_RATE)
        self.assertAlmostEqual(
            result["totGlycanConc"],
            api.DEFAULT_PROTEIN_PROD_RATE * api.DEFAULT_TAU / api.DEFAULT_COMPARTMENT_VOLUME,
        )

    def test_predict_applies_scalar_overrides(self):
        with mock.patch.object(api, "_run_model", return_value={"total": 1.0, "top": []}) as solver:
            result = api.predict(
                self.baseline,
                tau=7.5,
                compartment_volume=3.0,
                protein_prod_rate=500.0,
            )
        _, _, _, tau_v, vol_v, rate_v, _ = solver.call_args.args
        self.assertEqual(tau_v, 7.5)
        self.assertEqual(vol_v, 3.0)
        self.assertEqual(rate_v, 500.0)
        self.assertAlmostEqual(result["totGlycanConc"], 500.0 * 7.5 / 3.0)

    def test_predict_rejects_out_of_range_scalars(self):
        bad_calls = [
            ({"tau": -1.0},),
            ({"tau": 1000.0},),
            ({"compartment_volume": 0.0},),
            ({"protein_prod_rate": float("nan")},),
            ({"protein_prod_rate": "fast"},),
        ]
        for kwargs in bad_calls:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    api.predict(self.baseline, **kwargs[0])

    def test_predict_validates_donor_concs(self):
        with mock.patch.object(api, "_run_model", return_value={"total": 1.0, "top": []}):
            result = api.predict(self.baseline, donor_concs={"UDP-Gal": 1234.0})
        self.assertEqual(result["donorConcs"]["UDP-Gal"], 1234.0)
        self.assertEqual(
            result["donorConcs"]["UDP-GlcNAc"], api.DEFAULT_DONOR_CONC["UDP-GlcNAc"]
        )

        for bad in (-1, 0, float("inf"), True, "9200"):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    api.predict(self.baseline, donor_concs={"UDP-Gal": bad})

        with self.assertRaises(ValueError):
            api.predict(self.baseline, donor_concs={"H20": 5.0})

    def test_identical_predictions_use_bounded_cache(self):
        fake_result = {"total": 1.0, "top": []}
        with mock.patch.object(api, "_run_model", return_value=fake_result) as solver:
            first = api.predict(self.baseline)
            second = api.predict(self.baseline)
        self.assertFalse(first["cacheHit"])
        self.assertTrue(second["cacheHit"])
        solver.assert_called_once()

    def test_different_physiology_bypasses_cache(self):
        fake_result = {"total": 1.0, "top": []}
        with mock.patch.object(api, "_run_model", return_value=fake_result) as solver:
            api.predict(self.baseline)
            api.predict(self.baseline, tau=10.0)
            api.predict(self.baseline, donor_concs={"GDP-Fuc": 4000.0})
        self.assertEqual(solver.call_count, 3)


class StructureLabelTests(unittest.TestCase):
    def test_structure_file_covers_every_structure(self):
        tuples = api._read_structure_tuples(api._NETWORK_STRUCTURES_PATH)
        self.assertEqual(len(tuples), len(api._network.structures))
        for idx in api._network.structures:
            self.assertIn(int(idx), tuples)

    def test_file_tuples_match_in_memory_glycoforms(self):
        tuples = api._read_structure_tuples(api._NETWORK_STRUCTURES_PATH)
        for idx, glycoform in api._network.structures.items():
            self.assertEqual(tuple(glycoform), tuples[int(idx)])

    def test_build_structures_labels_come_from_file(self):
        from iupac import krambeck_to_iupac
        tuples = api._read_structure_tuples(api._NETWORK_STRUCTURES_PATH)
        records = {s["id"]: s for s in api._build_structures()}
        self.assertEqual(set(records), set(tuples))
        for struct_id, tup in tuples.items():
            self.assertEqual(records[struct_id]["label"], krambeck_to_iupac(tup))
            for field, value in zip(api._STRUCTURE_FIELDS, tup):
                self.assertEqual(records[struct_id][field], value)

    def test_high_mannose_reference_names(self):
        records = {s["id"]: s for s in api._build_structures()}
        self.assertEqual(
            records[0]["label"],
            "Mana1-2Mana1-6(Mana1-2Mana1-3)Mana1-6"
            "(Mana1-2Mana1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb1",
        )


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
        self.assertIn("donorDefaults", data)
        self.assertIn("tauDefault", data)
        self.assertIn("compartmentVolumeDefault", data)
        self.assertIn("proteinProdRateDefault", data)
        self.assertIn("totGlycanConcDefault", data)
        self.assertNotIn("promoterLevels", data)
        self.assertNotIn("enzymePresets", data)

    def test_frontend_serves_four_compartment_controls(self):
        with self.client.get("/") as response:
            self.assertEqual(response.status_code, 200)
            page = response.get_data(as_text=True)
        self.assertIn("Enzyme compartment distribution", page)
        self.assertIn("config.compartments", page)
        self.assertNotIn("Promoter strength", page)
        self.assertNotIn("Enzyme distribution presets", page)

    def test_missing_distribution_is_rejected_before_prediction(self):
        with mock.patch.object(flask_app, "predict") as predictor:
            response = self.client.post(
                "/api/predict",
                json={"tau": 5.56},
            )
        self.assertEqual(response.status_code, 400)
        self.assertIn("enzymeDistribution is required", response.get_json()["error"])
        predictor.assert_not_called()

    def test_invalid_distribution_is_rejected_before_solver(self):
        with mock.patch.object(api, "_run_model") as solver:
            response = self.client.post(
                "/api/predict",
                json={"enzymeDistribution": {"ManI": {}}},
            )
        self.assertEqual(response.status_code, 400)
        solver.assert_not_called()

    def test_invalid_physiology_is_rejected_before_solver(self):
        distribution = api._matrix_to_distribution(api._BASELINE_DIST_MATRIX)
        with mock.patch.object(api, "_run_model") as solver:
            response = self.client.post(
                "/api/predict",
                json={"enzymeDistribution": distribution, "tau": -1.0},
            )
        self.assertEqual(response.status_code, 400)
        solver.assert_not_called()

    def test_valid_prediction_returns_applied_parameters(self):
        distribution = api._matrix_to_distribution(api._BASELINE_DIST_MATRIX * 100)
        with api._PREDICTION_LOCK:
            api._PREDICTION_CACHE.clear()
        with mock.patch.object(api, "_run_model", return_value={"total": 1.0, "top": []}):
            response = self.client.post(
                "/api/predict",
                json={
                    "enzymeDistribution": distribution,
                    "tau": 7.0,
                    "proteinProdRate": 800.0,
                },
            )
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertIn("enzymeDistribution", data)
        self.assertFalse(data["cacheHit"])
        self.assertEqual(data["tau"], 7.0)
        self.assertEqual(data["proteinProdRate"], 800.0)
        self.assertAlmostEqual(
            data["totGlycanConc"], 800.0 * 7.0 / api.DEFAULT_COMPARTMENT_VOLUME
        )
        for profile in data["enzymeDistribution"].values():
            self.assertAlmostEqual(sum(profile.values()), 1)


if __name__ == "__main__":
    unittest.main()
