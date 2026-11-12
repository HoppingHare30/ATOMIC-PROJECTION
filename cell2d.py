"""
2D primitive in-plane lattice vectors, Gauss/Lagrange reduction, and cell area.

Identifies the 2D primitive unit cell of a crystallographic plane using exact
Bravais lattice parity rules, applies Minkowski-Gauss reduction, and computes cell areas.
"""
from __future__ import annotations

from fractions import Fraction
import math

from lattice_generator.core.plane import validate_and_reduce


def _dot_frac(
    u: tuple[Fraction, Fraction, Fraction],
    v: tuple[Fraction, Fraction, Fraction],
) -> Fraction:
    """Exact dot product of two 3-tuples of Fractions."""
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _cross_frac(
    u: tuple[Fraction, Fraction, Fraction],
    v: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    """Exact 3D cross product of two 3-tuples of Fractions."""
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def gauss_reduce(
    v1: tuple[Fraction, Fraction, Fraction],
    v2: tuple[Fraction, Fraction, Fraction],
) -> tuple[tuple[Fraction, Fraction, Fraction], tuple[Fraction, Fraction, Fraction]]:
    """
    Gauss/Lagrange (Minkowski) reduction of a 2D lattice basis.

    Implements: MATH_SPEC §7.2:
      While |v2|^2 < |v1|^2: swap(v1, v2)
      n = ceil(v1.v2 / |v1|^2 - 1/2) so that (v1.v2)/|v1|^2 in (-1/2, 1/2]
      v2 = v2 - n*v1
    Units: Vectors in units of lattice parameter a.
    Exact: Exact rational arithmetic via fractions.Fraction.

    Args:
        v1: First basis vector as Fraction 3-tuple.
        v2: Second basis vector as Fraction 3-tuple.

    Returns:
        Pair (v1_red, v2_red) of Minkowski-reduced Fraction 3-tuples satisfying
        |v1|² <= |v2|² and (v1.v2)/|v1|² in (-0.5, 0.5].
    """
    while True:
        if _dot_frac(v2, v2) < _dot_frac(v1, v1):
            v1, v2 = v2, v1
        len_sq1 = _dot_frac(v1, v1)
        if len_sq1 == Fraction(0, 1):
            raise ValueError("Zero vector encountered in Gauss reduction.")
        ratio = _dot_frac(v1, v2) / len_sq1
        # Shift ratio into (-1/2, 1/2]:
        n = math.ceil(ratio - Fraction(1, 2))
        if n != 0:
            fn = Fraction(n, 1)
            v2 = (v2[0] - fn * v1[0], v2[1] - fn * v1[1], v2[2] - fn * v1[2])

        if _dot_frac(v2, v2) < _dot_frac(v1, v1):
            v1, v2 = v2, v1
            continue
        break

    return (v1, v2)


def find_inplane_vectors(
    structure: str,
    h: int,
    k: int,
    l: int,
    bound: int = 4,
) -> tuple[tuple[Fraction, Fraction, Fraction], tuple[Fraction, Fraction, Fraction]]:
    """
    Find the two shortest independent in-plane primitive Bravais lattice vectors.

    Implements: MATH_SPEC §7.1 candidate enumeration with doubled-coord parity filters:
      SC:  all w_i even
      FCC: (w1 + w2 + w3) % 2 == 0 [sum-even]
      BCC: w1 % 2 == w2 % 2 == w3 % 2 [all same parity]
      hw1 + kw2 + lw3 == 0 (perpendicularity)
    followed by Gauss/Lagrange reduction (§7.2).
    Units: In-plane lattice vectors in units of lattice parameter a.
    Exact: Exact integer enumeration and Fraction arithmetic.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.
        bound: Initial search half-width in conventional cell units (default: 4).

    Returns:
        Pair (v1, v2) of Minkowski-reduced Fraction 3-tuples spanning the 2D primitive cell.

    Raises:
        ValueError: If structure is unknown, plane is invalid, or no independent pair is found.
    """
    h_red, k_red, l_red = validate_and_reduce(h, k, l)
    s = structure.strip().lower()

    if s not in ("fcc", "bcc", "sc"):
        raise ValueError(f"Unknown structure '{structure}'. Supported: 'fcc', 'bcc', 'sc'.")

    # Search candidates with expanding radius if necessary
    current_bound = bound
    while current_bound <= max(bound, 12):
        candidates: list[tuple[int, int, int]] = []
        max_w = 2 * current_bound
        for w1 in range(-max_w, max_w + 1):
            for w2 in range(-max_w, max_w + 1):
                for w3 in range(-max_w, max_w + 1):
                    if w1 == 0 and w2 == 0 and w3 == 0:
                        continue
                    # Parity filter per §7.1
                    if s == "sc":
                        if w1 % 2 != 0 or w2 % 2 != 0 or w3 % 2 != 0:
                            continue
                    elif s == "fcc":
                        if (w1 + w2 + w3) % 2 != 0:
                            continue
                    elif s == "bcc":
                        if not (w1 % 2 == w2 % 2 == w3 % 2):
                            continue

                    # Perpendicularity condition hw1 + kw2 + lw3 == 0
                    if h_red * w1 + k_red * w2 + l_red * w3 == 0:
                        candidates.append((w1, w2, w3))

        if len(candidates) >= 2:
            # Sort by squared length: w1^2 + w2^2 + w3^2
            candidates.sort(key=lambda w: w[0] * w[0] + w[1] * w[1] + w[2] * w[2])

            v1_cand = (
                Fraction(candidates[0][0], 2),
                Fraction(candidates[0][1], 2),
                Fraction(candidates[0][2], 2),
            )

            # Find shortest independent vector v2
            v2_cand: tuple[Fraction, Fraction, Fraction] | None = None
            for w in candidates[1:]:
                test_v2 = (Fraction(w[0], 2), Fraction(w[1], 2), Fraction(w[2], 2))
                cross = _cross_frac(v1_cand, test_v2)
                if cross != (Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)):
                    v2_cand = test_v2
                    break

            if v2_cand is not None:
                return gauss_reduce(v1_cand, v2_cand)

        current_bound += 2

    raise ValueError(f"Could not find two independent in-plane vectors for {structure} ({h},{k},{l}) within bound {current_bound}.")


def cell_area_sq(
    v1: tuple[Fraction, Fraction, Fraction],
    v2: tuple[Fraction, Fraction, Fraction],
) -> Fraction:
    """
    Compute exact rational square of 2D primitive cell area: A^2 = |v1 x v2|^2.

    Implements: MATH_SPEC §7.4: A_sq = (v1.v1)(v2.v2) - (v1.v2)^2.
    Units: Units of a^4.
    Exact: Exact rational arithmetic via fractions.Fraction.

    Args:
        v1: First in-plane vector as Fraction 3-tuple.
        v2: Second in-plane vector as Fraction 3-tuple.

    Returns:
        Fraction representing the exact squared area A^2.
    """
    g11 = _dot_frac(v1, v1)
    g12 = _dot_frac(v1, v2)
    g22 = _dot_frac(v2, v2)
    return g11 * g22 - g12 * g12


def cell_area(
    v1: tuple[Fraction, Fraction, Fraction],
    v2: tuple[Fraction, Fraction, Fraction],
) -> float:
    """
    Compute 2D primitive cell area |v1 x v2| as a float.

    Implements: MATH_SPEC §7.3: A_2D = |v1 x v2|.
    Units: Units of a^2.
    Exact: Square root of exact cell_area_sq.

    Args:
        v1: First in-plane vector as Fraction 3-tuple.
        v2: Second in-plane vector as Fraction 3-tuple.

    Returns:
        Float area of 2D primitive unit cell.
    """
    return math.sqrt(float(cell_area_sq(v1, v2)))
