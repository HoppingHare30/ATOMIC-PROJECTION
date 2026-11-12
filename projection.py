"""
Orthogonal projection onto crystallographic planes.

Calculates exact signed distance quantities, constructs right-handed
orthonormal in-plane frames (e1, e2, n_hat), and projects 3D atomic coordinates.
"""
from __future__ import annotations

from fractions import Fraction
import numpy as np

from lattice_generator.core.plane import unit_normal, validate_and_reduce


def signed_distance_exact(
    pos_frac: tuple[Fraction, Fraction, Fraction],
    h: int,
    k: int,
    l: int,
) -> Fraction:
    """
    Compute exact un-normalized signed distance h*x + k*y + l*z as a Fraction.

    Implements: MATH_SPEC §4: d' = hx + ky + lz.
    Units: Units of a (proportional to physical signed distance d * sqrt(h^2+k^2+l^2)).
    Exact: Exact rational arithmetic via fractions.Fraction.

    Args:
        pos_frac: Atomic position as 3-tuple of Fractions in units of a.
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.

    Returns:
        Fraction equal to h*pos_frac[0] + k*pos_frac[1] + l*pos_frac[2].
    """
    h_red, k_red, l_red = validate_and_reduce(h, k, l)
    return Fraction(h_red) * pos_frac[0] + Fraction(k_red) * pos_frac[1] + Fraction(l_red) * pos_frac[2]


def orthonormal_basis(
    h: int,
    k: int,
    l: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Construct right-handed orthonormal in-plane basis vectors (e1, e2, n_hat).

    Implements: MATH_SPEC §6.1 Gram-Schmidt construction using the smallest component of n_hat.
    Units: Dimensionless orthonormal unit vectors.
    Exact: Float64 vectors with orthonormality |e1|=|e2|=1, e1.e2=0, e1 x e2 = n_hat to 1e-12.

    Args:
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.

    Returns:
        Tuple (e1, e2, n_hat) of NumPy arrays of shape (3,).

    Raises:
        ValueError: If (h, k, l) == (0, 0, 0).
    """
    n_hat = unit_normal(h, k, l)

    # Choose non-parallel trial vector t along axis of smallest component
    ax, ay, az = abs(n_hat[0]), abs(n_hat[1]), abs(n_hat[2])
    if ax <= ay and ax <= az:
        t = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    elif ay <= az:
        t = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    else:
        t = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    # Gram-Schmidt to find e1
    e1_raw = t - np.dot(t, n_hat) * n_hat
    e1 = e1_raw / np.linalg.norm(e1_raw)

    # Right-hand rule: e2 = n_hat x e1 so that e1 x e2 = n_hat
    e2 = np.cross(n_hat, e1)
    e2 = e2 / np.linalg.norm(e2)

    return (e1, e2, n_hat)


def project_2d(
    positions: np.ndarray,
    e1: np.ndarray,
    e2: np.ndarray,
) -> np.ndarray:
    """
    Project 3D Cartesian coordinates onto the 2D plane defined by basis (e1, e2).

    Implements: MATH_SPEC §6.2: r_perp = (r.e1, r.e2).
    Units: Input and output positions in units of lattice parameter a.
    Exact: Standard float64 dot product.

    Args:
        positions: NumPy array of shape (N, 3) or (3,).
        e1: First in-plane unit basis vector of shape (3,).
        e2: Second in-plane unit basis vector of shape (3,).

    Returns:
        NumPy array of shape (N, 2) or (2,) with projected 2D coordinates.
    """
    pos = np.asarray(positions, dtype=np.float64)
    if pos.ndim == 1:
        return np.array([np.dot(pos, e1), np.dot(pos, e2)], dtype=np.float64)
    proj_x = np.dot(pos, e1)
    proj_y = np.dot(pos, e2)
    return np.column_stack([proj_x, proj_y])
