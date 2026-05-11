import math
import tempfile
import unittest
from pathlib import Path

import mitsuba as mi
import numpy as np

mi.set_variant("scalar_rgb")

import RealisticCamera  # noqa: E402
from utils import intersect_spherical_element, load_lens_file, refract  # noqa: E402


LENS = "scenes/kitchen/lenses/fisheye.10mm.dat"


def make_sensor(**kwargs):
    props = {
        "type": "scene",
        "sensor": {
            "type": "realisticcamera",
            "lens_file": LENS,
            "aperture_diameter": 6.0,
            "focus_distance": 3.0,
            "film_diagonal": 35.0,
            "mm_to_world": 0.001,
            "exit_pupil_sample_count": 256,
            "sampler": {"type": "independent", "sample_count": 4},
            "film": {"type": "hdrfilm", "width": 64, "height": 36},
        },
    }
    props["sensor"].update(kwargs)
    return mi.load_dict(props).sensors()[0]


def finite_vector(v):
    return all(math.isfinite(float(v[i])) for i in range(3))


class TestRealisticCamera(unittest.TestCase):
    def test_lens_file_loader_and_validation(self):
        values = load_lens_file(LENS)
        self.assertEqual(len(values), 48)
        np.testing.assert_allclose(values[:4], [30.2249, 0.8335, 1.62, 30.34])

        with tempfile.TemporaryDirectory() as tmpdir:
            bad = Path(tmpdir) / "bad.dat"
            bad.write_text("1 2 3\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_lens_file(bad)

    def test_aperture_clamping_and_film_extent(self):
        sensor = make_sensor(aperture_diameter=8.0)
        stop = next(e for e in sensor.elements if e.curvature_radius == 0)
        self.assertAlmostEqual(stop.aperture_radius, 0.5 * 6.08e-3)
        expected_width = 0.035 / math.sqrt(1.0 + (36 / 64) ** 2)
        self.assertAlmostEqual(
            sensor.physical_extent.max_x - sensor.physical_extent.min_x,
            expected_width,
        )

    def test_intersection_refraction_and_tir(self):
        sensor = make_sensor()
        hit = intersect_spherical_element(1.0, 2.0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
        self.assertIsNotNone(hit)
        self.assertAlmostEqual(hit[0], 1.0)

        wt = refract((0.0, 0.0, -1.0), (0.0, 0.0, 1.0), 1.5)
        np.testing.assert_allclose(wt, (0.0, 0.0, 1.0))
        self.assertIsNone(
            refract((math.sqrt(0.75), 0.0, 0.5), (0.0, 0.0, 1.0), 0.5)
        )

    def test_xml_load_and_sample_rays_are_finite(self):
        sensor = make_sensor()
        center_ray, center_w = sensor.sample_ray(
            0.0, 0.5, mi.Point2f(32.5, 18.5), mi.Point2f(0.5, 0.5)
        )
        edge_ray, _ = sensor.sample_ray(0.0, 0.5, mi.Point2f(8.0, 8.0), mi.Point2f(0.5, 0.5))

        self.assertGreater(center_w[0], 0.0)
        self.assertTrue(finite_vector(center_ray.o))
        self.assertTrue(finite_vector(center_ray.d))
        self.assertTrue(finite_vector(edge_ray.o))
        self.assertTrue(finite_vector(edge_ray.d))

    def test_low_res_render_smoke(self):
        scene = mi.load_dict(
            {
                "type": "scene",
                "integrator": {"type": "path"},
                "sensor": {
                    "type": "realisticcamera",
                    "lens_file": LENS,
                    "aperture_diameter": 6.0,
                    "focus_distance": 3.0,
                    "film_diagonal": 35.0,
                    "mm_to_world": 0.001,
                    "exit_pupil_sample_count": 128,
                    "sampler": {"type": "independent", "sample_count": 2},
                    "film": {"type": "hdrfilm", "width": 16, "height": 9},
                },
                "env": {
                    "type": "constant",
                    "radiance": {"type": "rgb", "value": [1.0, 1.0, 1.0]},
                },
            }
        )
        try:
            image = mi.render(scene, spp=2)
        except RuntimeError as exc:
            if "LLVM-C.dll" in str(exc):
                self.skipTest("scalar Mitsuba render backend is unavailable")
            raise
        self.assertGreater(float(np.array(image).mean()), 0.0)


if __name__ == "__main__":
    unittest.main()
