"""S-JTSK / Krovak East North (sjtsk.py) against PROJ 9 (pyproj 3.8, EPSG:5514 with towgs84
570.8,85.7,462.8,4.998,1.587,5.261,3.56; the points below computed with it, 2026-10-30)."""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import sjtsk  # noqa: E402

# lat, lon, x, y: Prague (Nusle), Brno, Liberec, the Austrian border, Ostrava
PROJ = [
    (50.0645, 14.4296, -742565.931, -1045624.640),
    (49.1951, 16.6068, -598248.812, -1160744.686),
    (50.7663, 15.0543, -688246.909, -974109.515),
    (48.8127, 14.3175, -769682.783, -1182415.697),
    (49.8209, 18.2625, -472181.983, -1103062.861),
]


class SjtskTest(unittest.TestCase):
    def test_forward(self):
        la, lo, x, y = map(np.array, zip(*PROJ))
        fx, fy = sjtsk.forward(la, lo)
        np.testing.assert_allclose(fx, x, atol=0.01)
        np.testing.assert_allclose(fy, y, atol=0.01)

    def test_inverse(self):
        la, lo, x, y = map(np.array, zip(*PROJ))
        ila, ilo = sjtsk.inverse(x, y)
        np.testing.assert_allclose(ila, la, atol=1e-7)     # (~1 cm)
        np.testing.assert_allclose(ilo, lo, atol=1e-7)

    def test_scalar(self):
        x, y = sjtsk.forward(50.0645, 14.4296)
        self.assertAlmostEqual(float(x), PROJ[0][2], delta=0.01)
        self.assertAlmostEqual(float(y), PROJ[0][3], delta=0.01)


if __name__ == "__main__":
    unittest.main()
