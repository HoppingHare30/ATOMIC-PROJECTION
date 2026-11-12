"""
Interplanar spacing and planar density computations.

Computes theoretical formula spacing d_formula, true atomic interplanar spacing
d_true measured from atomic layer positions (accounting for extinction), and
planar density rho = N_conv * d_true.
"""
from __future__ import annotations

from fractions import Fraction
import math

from lattice_generator.core.lattice import get_basis, generate_supercell_exact


def formula_spacing(h: int, k: int, l: int) -> float:
    """
    Compute theoretical interplanar spacing from standard crystallographic formula.

    Implements: MATH_SPEC §5.1: d_formula = a / sqrt(h^2 + k^2 + l^2) with a=1.
    Units: Spacing in units of lattice parameter a.
    Exact: Computed using float64 square root.

    Args:
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.

    Returns:
        Float interplanar spacing d_formula in units of a.

    Raises:
        ValueError: If (h, k, l) == (0, 0, 0).
    """
    if h == 0 and k == 0 and l == 0:
        raise ValueError("Miller indices (h, k, l) cannot all be zero.")
    return 1.0 / math.sqrt(h * h + k * k + l * l)


def _get_nconv(structure: str) -> int:
    """Return number of atoms per conventional cubic cell."""
    s = structure.strip().lower()
    if s == "fcc":
        return 4
    elif s == "bcc":
        return 2
    elif s == "sc":
        return 1
    else:
        raise ValueError(f"Unknown structure '{structure}'. Supported: 'fcc', 'bcc', 'sc'.")


def true_spacing_sq(
    structure: str,
    h: int,
    k: int,
    l: int,
    N: int = 2,
) -> Fraction:
    """
    Compute exact rational square of true atomic interplanar spacing (d_true)^2.

    Implements: MATH_SPEC §5.2–5.3, §7.4:
      (d_true)^2 = (Delta_d_int)^2 / (4 * (h^2 + k^2 + l^2))
    where Delta_d_int is the minimum integer difference between distinct 2*(hx+ky+lz) layer keys.
    Units: Units of a^2.
    Exact: Exact rational arithmetic via fractions.Fraction.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.
        N: Supercell radius for measuring layer spacing (default: 2 is sufficient).

    Returns:
        Fraction representing the exact squared true spacing (d_true)^2.

    Raises:
        ValueError: If (h, k, l) == (0, 0, 0) or structure is unknown.
    """
    if h == 0 and k == 0 and l == 0:
        raise ValueError("Miller indices (h, k, l) cannot all be zero.")

    # Generate small centered supercell
    radius = max(N, 2)
    atoms = generate_supercell_exact(structure, radius)

    # Compute 2*(hx + ky + lz) as exact integers
    fh = Fraction(h, 1)
    fk = Fraction(k, 1)
    fl = Fraction(l, 1)

    layer_keys_set: set[int] = set()
    for pos in atoms:
        val = Fraction(2, 1) * (fh * pos[0] + fk * pos[1] + fl * pos[2])
        layer_keys_set.add(int(val))

    sorted_keys = sorted(layer_keys_set)
    if len(sorted_keys) < 2:
        raise ValueError("Supercell too small to determine layer spacing.")

    diffs = [sorted_keys[i + 1] - sorted_keys[i] for i in range(len(sorted_keys) - 1)]
    delta_d_int = min(diffs)

    h_sq_sum = h * h + k * k + l * l
    return Fraction(delta_d_int * delta_d_int, 4 * h_sq_sum)


def true_spacing(
    structure: str,
    h: int,
    k: int,
    l: int,
    N: int = 2,
) -> float:
    """
    Compute true atomic interplanar spacing d_true from atomic layer positions.

    Implements: MATH_SPEC §5.2–5.3.
    Units: Spacing in units of lattice parameter a.
    Exact: Square root of exact rational true_spacing_sq.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.
        N: Supercell radius for measuring layer spacing.

    Returns:
        Float true interplanar spacing d_true in units of a.
    """
    return math.sqrt(float(true_spacing_sq(structure, h, k, l, N)))


def planar_density_sq(
    structure: str,
    h: int,
    k: int,
    l: int,
    N: int = 2,
) -> Fraction:
    """
    Compute exact rational square of planar density: rho^2 = N_conv^2 * (d_true)^2.

    Implements: MATH_SPEC §8.1, §7.4.
    Units: Units of a^(-4).
    Exact: Exact rational arithmetic via fractions.Fraction.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.
        N: Supercell radius for measuring layer spacing.

    Returns:
        Fraction representing rho^2.
    """
    n_conv = _get_nconv(structure)
    d_sq = true_spacing_sq(structure, h, k, l, N)
    return Fraction(n_conv * n_conv, 1) * d_sq


def planar_density(
    structure: str,
    h: int,
    k: int,
    l: int,
    N: int = 2,
) -> float:
    """
    Compute planar atomic density rho = N_conv * d_true.

    Implements: MATH_SPEC §8.1: rho_hkl = N_conv * d_true / a^3 (with a=1).
    Units: Atoms per unit area in units of a^(-2).
    Exact: Evaluated from exact planar_density_sq.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.
        N: Supercell radius for measuring layer spacing.

    Returns:
        Float planar density in units of a^(-2).
    """
    return math.sqrt(float(planar_density_sq(structure, h, k, l, N)))
