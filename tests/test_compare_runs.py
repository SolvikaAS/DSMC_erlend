"""Checks for profile extraction and matching saved comparison cases."""

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "postprocessing"))
from compare_runs import check_settings, extract_profile, load_run, load_sparta, resolve_field, run_label


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        rows = [[x, y, 1e17, 10, 20, 200, 200, 200]
                for x in (-4, 0, 4) for y in (-4, 0, 4)]
        self.run = {
            "path": Path("case.plt"), "data": np.array(rows, dtype=float),
            "metadata": {
                "obstacle": {"type": "box", "center": [0, 0], "width": 2, "height": 2},
                "wall_model": "specular", "angle_of_attack_degrees": 0,
                "freestream": {"number_density": 1e17},
            },
        }

    def test_horizontal_profile_and_solid_gap(self):
        coordinates, values, actual = extract_profile(self.run, "x", 0.2)
        np.testing.assert_array_equal(coordinates, [-4, 0, 4])
        self.assertEqual(actual, 0)
        self.assertTrue(np.isnan(values["temperature"][1]))
        self.assertEqual(values["temperature"][0], 200)
        self.assertAlmostEqual(values["pressure"][0], 1e17 * 1.380662e-23 * 200)
        self.assertAlmostEqual(values["speed"][0], np.hypot(10, 20))

    def test_vertical_profile(self):
        coordinates, values, actual = extract_profile(self.run, "y", -3)
        self.assertEqual(actual, -4)
        np.testing.assert_array_equal(coordinates, [-4, 0, 4])
        self.assertTrue(np.isfinite(values["temperature"]).all())

    def test_wall_models_can_differ_but_density_must_match(self):
        diffuse = copy.deepcopy(self.run)
        diffuse["metadata"]["wall_model"] = "diffuse"
        check_settings(self.run, diffuse)
        diffuse["metadata"]["freestream"]["number_density"] = 1e14
        with self.assertRaisesRegex(ValueError, "freestream"):
            check_settings(self.run, diffuse)

    def test_empty_sampled_profile_is_rejected(self):
        self.run["data"][:, 2] = 0
        with self.assertRaisesRegex(ValueError, "No sampled gas cells"):
            extract_profile(self.run, "x", 0)

    def test_cll_case_can_be_compared_and_labelled(self):
        cll = copy.deepcopy(self.run)
        cll["metadata"]["wall_model"] = "cll"
        cll["metadata"]["cll_accommodation"] = {"normal": 0.8, "tangential": 0.6}
        check_settings(self.run, cll)
        self.assertEqual(run_label(cll), "CLL (an=0.8, at=0.6)")

    def test_saved_field_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            field = Path(directory) / "specular" / "DSMC2D_Field.plt"
            field.parent.mkdir()
            np.savetxt(field, self.run["data"], header="VARIABLES\nZONE", comments="")
            field.with_suffix(".json").write_text(json.dumps(self.run["metadata"]))
            self.assertEqual(resolve_field(field.parent), field)
            loaded = load_run(field.parent)
            np.testing.assert_array_equal(loaded["data"], self.run["data"])

    def test_sparta_profile_without_sidecar_preserves_density(self):
        with tempfile.TemporaryDirectory() as directory:
            field = Path(directory) / "sparta.plt"
            data = self.run["data"][self.run["data"][:, 0] == -4].copy()
            data[:, 2] = 1e18
            np.savetxt(field, data, header="VARIABLES\nZONE", comments="")
            run = load_sparta(field)
            coordinates, values, actual = extract_profile(run, "y", -4)
            np.testing.assert_array_equal(coordinates, [-4, 0, 4])
            np.testing.assert_array_equal(values["density"], [1e18, 1e18, 1e18])
            self.assertEqual(actual, -4)
            self.assertEqual(run["label"], "SPARTA (diffuse)")
            with self.assertRaisesRegex(ValueError, "update PROFILE_POSITION"):
                extract_profile(run, "y", -3)
            with self.assertRaisesRegex(ValueError, "no profile along x"):
                extract_profile(run, "x", 0)


if __name__ == "__main__":
    unittest.main()
