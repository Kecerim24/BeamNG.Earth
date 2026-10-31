"""The Czech Republic's S-JTSK / Krovak East North (EPSG:5514), pure numpy (no PROJ in the app): ČÚZK's elevation data
(roadz.py: DMP 1G for the Czech bridge decks) is in it.

    x, y = forward(lat, lon)            # WGS84 degrees (arrays) -> EPSG:5514 easting / northing (m, both negative)
    lat, lon = inverse(x, y)            # and back

The Krovak oblique conformal conic on Bessel 1841 (EPSG Guidance Note 7-2, method 1041), the datum by the 7-parameter
shift PROJ gives EPSG:5514 (towgs84 570.8,85.7,462.8,4.998,1.587,5.261,3.56, position vector; ~1 m, against ~10 m
for the old 3-parameter one). Heights are left alone (h = 0 for the datum shift: < 1 cm here). Checked against PROJ
in tests/test_sjtsk.py.
"""
import math

import numpy as np

A = 6377397.155                     # Bessel 1841
F = 1 / 299.1528128
E2 = F * (2 - F)
E = math.sqrt(E2)
A84 = 6378137.0
F84 = 1 / 298.257223563
E2_84 = F84 * (2 - F84)

PHI_C = math.radians(49.5)                                  # latitude of the projection centre
LAM_0 = math.radians(24 + 50 / 60)                          # its longitude (42 30' east of Ferro)
ALPHA_C = math.radians(30 + 17 / 60 + 17.30311 / 3600)      # co-latitude of the cone axis
PHI_P = math.radians(78.5)                                  # the pseudo standard parallel
K_P = 0.9999

TOWGS84 = (570.8, 85.7, 462.8, 4.998, 1.587, 5.261, 3.56)   # m, arc-seconds, ppm

_A = A * math.sqrt(1 - E2) / (1 - E2 * math.sin(PHI_C) ** 2)
_B = math.sqrt(1 + E2 * math.cos(PHI_C) ** 4 / (1 - E2))
_G0 = math.asin(math.sin(PHI_C) / _B)
_T0 = (math.tan(math.pi / 4 + _G0 / 2) * ((1 + E * math.sin(PHI_C)) / (1 - E * math.sin(PHI_C))) ** (E * _B / 2)
       / math.tan(math.pi / 4 + PHI_C / 2) ** _B)
_N = math.sin(PHI_P)
_R0 = K_P * _A / math.tan(PHI_P)


def _rot():
    rx, ry, rz = (math.radians(v / 3600) for v in TOWGS84[3:6])
    return np.array([[1.0, -rz, ry], [rz, 1.0, -rx], [-ry, rx, 1.0]])


_R = _rot()
_S = 1 + TOWGS84[6] * 1e-6
_T = np.array(TOWGS84[:3])


def _ecef(phi, lam, a, e2):
    n = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)
    return np.stack([n * np.cos(phi) * np.cos(lam), n * np.cos(phi) * np.sin(lam), n * (1 - e2) * np.sin(phi)], -1)


def _geodetic(xyz, a, e2):
    """(phi, lam) radians of ECEF points xyz (..., 3) (Bowring's, then two refinements: sub-mm near the ground)."""
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    p = np.hypot(x, y)
    phi = np.arctan2(z, p * (1 - e2))
    for _ in range(4):
        n = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)
        h = p / np.cos(phi) - n
        phi = np.arctan2(z, p * (1 - e2 * n / (n + h)))
    return phi, np.arctan2(y, x)


def wgs84_to_bessel(lat, lon):
    """(phi, lam) radians on Bessel / S-JTSK of WGS84 degrees."""
    xyz = _ecef(np.radians(np.asarray(lat, np.float64)), np.radians(np.asarray(lon, np.float64)), A84, E2_84)
    return _geodetic(((xyz - _T) / _S) @ _R, A, E2)        # (R's inverse: its transpose, the angles being small)


def bessel_to_wgs84(phi, lam):
    """WGS84 (lat, lon) degrees of Bessel / S-JTSK (phi, lam) radians."""
    xyz = _T + _S * (_ecef(phi, lam, A, E2) @ _R.T)
    p, l_ = _geodetic(xyz, A84, E2_84)
    return np.degrees(p), np.degrees(l_)


def forward(lat, lon):
    """EPSG:5514 (x, y) m of WGS84 (lat, lon) degrees."""
    phi, lam = wgs84_to_bessel(lat, lon)
    es = E * np.sin(phi)
    u = 2 * (np.arctan(_T0 * np.tan(phi / 2 + np.pi / 4) ** _B / ((1 + es) / (1 - es)) ** (E * _B / 2)) - np.pi / 4)
    v = _B * (LAM_0 - lam)
    t = np.arcsin(math.cos(ALPHA_C) * np.sin(u) + math.sin(ALPHA_C) * np.cos(u) * np.cos(v))
    d = np.arcsin(np.cos(u) * np.sin(v) / np.cos(t))
    r = _R0 * math.tan(math.pi / 4 + PHI_P / 2) ** _N / np.tan(t / 2 + np.pi / 4) ** _N
    th = _N * d
    return -r * np.sin(th), -r * np.cos(th)                 # (Krovak's westing / southing, negated: East North)


def inverse(x, y):
    """WGS84 (lat, lon) degrees of EPSG:5514 (x, y) m."""
    xs, ys = -np.asarray(y, np.float64), -np.asarray(x, np.float64)     # southing, westing
    r = np.hypot(xs, ys)
    d = np.arctan2(ys, xs) / _N
    t = 2 * (np.arctan((_R0 / r) ** (1 / _N) * math.tan(math.pi / 4 + PHI_P / 2)) - np.pi / 4)
    u = np.arcsin(math.cos(ALPHA_C) * np.sin(t) - math.sin(ALPHA_C) * np.cos(t) * np.cos(d))
    v = np.arcsin(np.cos(t) * np.sin(d) / np.cos(u))
    lam = LAM_0 - v / _B
    k = _T0 ** (-1 / _B) * np.tan(u / 2 + np.pi / 4) ** (1 / _B)
    phi = u
    for _ in range(6):
        es = E * np.sin(phi)
        phi = 2 * (np.arctan(k * ((1 + es) / (1 - es)) ** (E / 2)) - np.pi / 4)
    return bessel_to_wgs84(phi, lam)
