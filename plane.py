"""
Miller indices validation, reduction, and unit normal calculation.

Provides exact integer reduction of (hkl) planes and unit normal vectors
for cubic crystal geometries.
"""
from __future__ import annotations

import math
import numpy as np


def validate_and_reduce(h: int, k: int, l: int) -> tuple[int, int, int]:
    """
    Validate Miller indices (h, k, l) and reduce by their greatest common divisor.

    Implements: MATH_SPEC §3.3: (h, k, l) reduced by g = gcd(|h|, |k|, |l|).
    Units: Dimensionless reciprocal lattice plane indices.
    Exact: Exact integer arithmetic via math.gcd.

    Args:
        h: Miller index h (integer).
        k: Miller index k (integer).
        l: Miller index l (integer).

    Returns:
        3-tuple (h_red, k_red, l_red) of reduced integers.

    Raises:
        ValueError: If (h, k, l) == (0, 0, 0).
    """
    if h == 0 and k == 0 and l == 0:
        raise ValueError("Miller indices (h, k, l) cannot all be zero; plane is undefined.")
    g = math.gcd(math.gcd(abs(h), abs(k)), abs(l))
    return (h // g, k // g, l // g)


def unit_normal(h: int, k: int, l: int) -> np.ndarray:
    """
    Compute unit normal vector n_hat = (h, k, l) / sqrt(h^2 + k^2 + l^2).

    Implements: MATH_SPEC §3.2: n_hat = G_{hkl} / |G_{hkl}|.
    Units: Dimensionless unit vector in Cartesian cubic basis.
    Exact: Uses float64 square root and normalization.

    Args:
        h: Miller index h (integer).
        k: Miller index k (integer).
        l: Miller index l (integer).

    Returns:
        NumPy array of shape (3,) with |n_hat| = 1.0.

    Raises:
        ValueError: If (h, k, l) == (0, 0, 0).
    """
    h_red, k_red, l_red = validate_and_reduce(h, k, l)
    norm = math.sqrt(h_red * h_red + k_red * k_red + l_red * l_red)
    return np.array([h_red / norm, k_red / norm, l_red / norm], dtype=np.float64)
