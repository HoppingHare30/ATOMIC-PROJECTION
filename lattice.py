"""
Bravais lattice bases and supercell generation.

Provides basis vectors for FCC, BCC, and SC lattices in exact fractional
coordinates, and generates symmetric centered supercells.
"""
from __future__ import annotations

from fractions import Fraction
import numpy as np


def get_basis(structure: str) -> list[tuple[Fraction, Fraction, Fraction]]:
    """
    Return fractional-coordinate basis atoms for the conventional unit cell.

    Implements: MATH_SPEC §1.2.
    Units: Fractional coordinates of conventional cubic lattice parameter a.
    Exact: Exact arithmetic using fractions.Fraction to eliminate rounding errors.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', or 'sc', case-insensitive).

    Returns:
        List of 3-tuples (u, v, w) of Fraction objects representing atomic coordinates.

    Raises:
        ValueError: If structure is not recognized.
    """
    s = structure.strip().lower()
    if s == "fcc":
        return [
            (Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
            (Fraction(1, 2), Fraction(1, 2), Fraction(0, 1)),
            (Fraction(1, 2), Fraction(0, 1), Fraction(1, 2)),
            (Fraction(0, 1), Fraction(1, 2), Fraction(1, 2)),
        ]
    elif s == "bcc":
        return [
            (Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
            (Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)),
        ]
    elif s == "sc":
        return [
            (Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
        ]
    else:
        raise ValueError(f"Unknown structure '{structure}'. Supported structures are 'fcc', 'bcc', 'sc'.")


def generate_supercell_exact(
    structure: str,
    N: int | tuple[int, int, int]
) -> list[tuple[Fraction, Fraction, Fraction]]:
    """
    Generate centered supercell atom positions in exact Fraction coordinates.

    Implements: MATH_SPEC §1.3: r_{(n1,n2,n3,b)} = (n1+tau_b^x, n2+tau_b^y, n3+tau_b^z)
    for n_i in {-N_i, ..., N_i}.
    Units: Position coordinates in units of lattice parameter a.
    Exact: Exact rational arithmetic via fractions.Fraction.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        N: Supercell radius or (N1, N2, N3) tuple.

    Returns:
        List of 3-tuples (x, y, z) in Fraction coordinates.

    Raises:
        ValueError: If structure is unknown or N contains negative bounds.
    """
    basis = get_basis(structure)
    if isinstance(N, int):
        if N < 0:
            raise ValueError(f"Supercell radius N must be non-negative, got {N}")
        n1_max = n2_max = n3_max = N
    else:
        n1_max, n2_max, n3_max = N
        if n1_max < 0 or n2_max < 0 or n3_max < 0:
            raise ValueError(f"Supercell bounds must be non-negative, got {N}")

    positions: list[tuple[Fraction, Fraction, Fraction]] = []
    for n1 in range(-n1_max, n1_max + 1):
        f_n1 = Fraction(n1, 1)
        for n2 in range(-n2_max, n2_max + 1):
            f_n2 = Fraction(n2, 1)
            for n3 in range(-n3_max, n3_max + 1):
                f_n3 = Fraction(n3, 1)
                for tau in basis:
                    positions.append((f_n1 + tau[0], f_n2 + tau[1], f_n3 + tau[2]))
    return positions


def generate_supercell(
    structure: str,
    N: int | tuple[int, int, int]
) -> np.ndarray:
    """
    Generate centered supercell atom positions as a float64 NumPy array.

    Implements: MATH_SPEC §1.3.
    Units: Position coordinates in units of lattice parameter a.
    Exact: Derived from exact Fraction coordinates, cast to float64 at boundary.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        N: Supercell radius or (N1, N2, N3) tuple.

    Returns:
        NumPy array of shape (M, 3) and dtype float64.

    Raises:
        ValueError: If structure is unknown or N is negative.
    """
    frac_positions = generate_supercell_exact(structure, N)
    coords = [[float(p[0]), float(p[1]), float(p[2])] for p in frac_positions]
    return np.array(coords, dtype=np.float64)
