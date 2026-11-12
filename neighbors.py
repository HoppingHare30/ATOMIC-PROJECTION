"""
Nearest-neighbour distance and coordination number calculations.

Computes d_NN and coordination number CN from geometry by building a centered
supercell and grouping distances using exact squared integer arithmetic.
"""
from __future__ import annotations

from fractions import Fraction
import math

from lattice_generator.core.lattice import get_basis


def nearest_neighbour(
    structure: str | list[tuple[Fraction, Fraction, Fraction]],
    N: int = 3,
) -> tuple[float, int]:
    """
    Compute nearest-neighbour distance d_NN and coordination number CN from geometry.

    Implements: MATH_SPEC §10.2:
      Build centered NxNxN supercell, compare 4*d^2 as exact integers to group
      nearest neighbours without hardcoding.
    Units: d_NN in units of lattice parameter a; CN is an exact dimensionless integer.
    Exact: Distances grouped by exact integer 4*d^2 = (2x)^2 + (2y)^2 + (2z)^2 in Z.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc') OR list of basis
                   coordinate 3-tuples (u, v, w) of Fractions.
        N: Supercell radius (default: 3 conventional cells).

    Returns:
        Tuple (d_NN, CN) where d_NN is a float in units of a and CN is an integer.

    Raises:
        ValueError: If structure is unknown or basis has no neighbours.
    """
    if isinstance(structure, str):
        basis = get_basis(structure)
    else:
        basis = structure

    if not basis:
        raise ValueError("Basis cannot be empty.")

    # Atoms in centered supercell
    # Origin atom is taken as basis atom 0 in cell (0, 0, 0)
    origin = basis[0]

    min_4d2: int | None = None
    cn_count = 0

    radius = max(N, 2)
    for n1 in range(-radius, radius + 1):
        f_n1 = Fraction(n1, 1)
        for n2 in range(-radius, radius + 1):
            f_n2 = Fraction(n2, 1)
            for n3 in range(-radius, radius + 1):
                f_n3 = Fraction(n3, 1)
                for tau in basis:
                    dx = f_n1 + tau[0] - origin[0]
                    dy = f_n2 + tau[1] - origin[1]
                    dz = f_n3 + tau[2] - origin[2]

                    if dx == Fraction(0, 1) and dy == Fraction(0, 1) and dz == Fraction(0, 1):
                        continue

                    # Integer 2*dx, 2*dy, 2*dz
                    int_2dx = int(Fraction(2, 1) * dx)
                    int_2dy = int(Fraction(2, 1) * dy)
                    int_2dz = int(Fraction(2, 1) * dz)

                    four_d2 = int_2dx * int_2dx + int_2dy * int_2dy + int_2dz * int_2dz

                    if min_4d2 is None or four_d2 < min_4d2:
                        min_4d2 = four_d2
                        cn_count = 1
                    elif four_d2 == min_4d2:
                        cn_count += 1

    if min_4d2 is None:
        raise ValueError("No neighbour atoms found in supercell.")

    d_nn = math.sqrt(min_4d2) / 2.0
    return (d_nn, cn_count)
