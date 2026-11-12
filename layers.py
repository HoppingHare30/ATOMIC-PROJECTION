"""
Atomic layer assignment, filtering, and stacking sequence identification.

Groups atoms into discrete crystallographic layers using exact integer keys,
filters slabs by layer count or physical distance, and identifies stacking sequences
(e.g., ABC for FCC (111), AB for BCC (110)) via exact 2D lattice modulo reduction.
"""
from __future__ import annotations

from fractions import Fraction
import math
import numpy as np

from lattice_generator.core.plane import validate_and_reduce


def _ensure_fraction_tuples(
    positions: list[tuple[Fraction, Fraction, Fraction]] | np.ndarray,
) -> list[tuple[Fraction, Fraction, Fraction]]:
    """Convert input positions to list of 3-tuples of Fractions."""
    if isinstance(positions, np.ndarray):
        result: list[tuple[Fraction, Fraction, Fraction]] = []
        for row in positions:
            result.append((
                Fraction(row[0]).limit_denominator(1000),
                Fraction(row[1]).limit_denominator(1000),
                Fraction(row[2]).limit_denominator(1000),
            ))
        return result
    elif isinstance(positions, list) and len(positions) > 0 and not isinstance(positions[0][0], Fraction):
        return [
            (
                Fraction(p[0]).limit_denominator(1000),
                Fraction(p[1]).limit_denominator(1000),
                Fraction(p[2]).limit_denominator(1000),
            )
            for p in positions
        ]
    return positions  # type: ignore[return-value]


def assign_layers(
    positions_frac: list[tuple[Fraction, Fraction, Fraction]] | np.ndarray,
    h: int,
    k: int,
    l: int,
) -> dict[int, list[int]]:
    """
    Group atom indices by exact layer key = 2*(h*x + k*y + l*z).

    Implements: MATH_SPEC §4: integer dot product 2d' = 2(hx + ky + lz) in Z.
    Units: Layer keys are exact dimensionless integers; positions in units of a.
    Exact: Exact integer arithmetic without float tolerances.

    Args:
        positions_frac: List of Fraction 3-tuples (or NumPy array) of atomic coordinates.
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.

    Returns:
        Dict mapping integer layer key -> list of atom indices belonging to that layer.
    """
    h_red, k_red, l_red = validate_and_reduce(h, k, l)
    frac_positions = _ensure_fraction_tuples(positions_frac)

    fh = Fraction(h_red, 1)
    fk = Fraction(k_red, 1)
    fl = Fraction(l_red, 1)

    layer_dict: dict[int, list[int]] = {}
    for idx, pos in enumerate(frac_positions):
        val = Fraction(2, 1) * (fh * pos[0] + fk * pos[1] + fl * pos[2])
        if val.denominator != 1:
            raise ValueError(f"Atom coordinate {pos} produced non-integer layer key: {val}")
        key = int(val)
        layer_dict.setdefault(key, []).append(idx)

    return layer_dict


def filter_by_n_layers(
    layer_dict: dict[int, list[int]],
    positions_frac: list[tuple[Fraction, Fraction, Fraction]] | np.ndarray,
    h: int,
    k: int,
    l: int,
    n_layers: int,
) -> dict[int, list[int]]:
    """
    Filter layer dictionary to retain at most n_layers on each side of reference plane (d=0).

    Implements: MATH_SPEC §1.3, §13.2; FEATURES F-05.
    Units: Layer counts (integers).
    Exact: Exact key selection based on sign of layer_key.

    Args:
        layer_dict: Dict mapping layer_key -> atom_indices.
        positions_frac: Atomic positions in units of a.
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.
        n_layers: Number of distinct layers to retain per side (d <= 0 and d > 0).

    Returns:
        Filtered dictionary mapping layer_key -> atom_indices.

    Raises:
        ValueError: If n_layers < 1 or if no layers match.
    """
    if n_layers < 1:
        raise ValueError(f"n_layers must be at least 1, got {n_layers}")

    sorted_keys = sorted(layer_dict.keys())
    neg_or_zero = [k for k in sorted_keys if k <= 0]
    pos = [k for k in sorted_keys if k > 0]

    selected_neg = neg_or_zero[-n_layers:] if len(neg_or_zero) > 0 else []
    selected_pos = pos[:n_layers] if len(pos) > 0 else []

    selected_keys = sorted(selected_neg + selected_pos)
    if not selected_keys:
        raise ValueError("No layers found for requested n_layers.")

    return {k: layer_dict[k] for k in selected_keys}


def filter_by_distance(
    layer_dict: dict[int, list[int]],
    positions_frac: list[tuple[Fraction, Fraction, Fraction]] | np.ndarray,
    h: int,
    k: int,
    l: int,
    d_min: float,
    d_max: float,
) -> dict[int, list[int]]:
    """
    Filter layers by signed physical distance d from reference plane: d_min <= d <= d_max.

    Implements: MATH_SPEC §4, §5: d = key / (2 * sqrt(h^2 + k^2 + l^2)).
    Units: Distances in units of lattice parameter a.
    Exact: Layer distance derived from exact integer key and plane norm.

    Args:
        layer_dict: Dict mapping layer_key -> atom_indices.
        positions_frac: Atomic positions in units of a.
        h: Miller index h.
        k: Miller index k.
        l: Miller index l.
        d_min: Minimum physical signed distance in units of a.
        d_max: Maximum physical signed distance in units of a.

    Returns:
        Filtered dictionary of layer_key -> atom_indices.

    Raises:
        ValueError: If d_min > d_max or no layers fall within the range.
    """
    if d_min > d_max:
        raise ValueError(f"d_min ({d_min}) cannot exceed d_max ({d_max})")

    h_red, k_red, l_red = validate_and_reduce(h, k, l)
    norm = math.sqrt(h_red * h_red + k_red * k_red + l_red * l_red)
    factor = 1.0 / (2.0 * norm)

    filtered: dict[int, list[int]] = {}
    for key, indices in layer_dict.items():
        d = key * factor
        if d_min - 1e-12 <= d <= d_max + 1e-12:
            filtered[key] = indices

    if not filtered:
        raise ValueError(f"No atoms found in distance range [{d_min}, {d_max}].")

    return filtered


def stacking_sequence(
    layers: dict[int, list[int]] | list[int],
    positions_frac: list[tuple[Fraction, Fraction, Fraction]] | np.ndarray,
    v1_frac: tuple[Fraction, Fraction, Fraction],
    v2_frac: tuple[Fraction, Fraction, Fraction],
) -> list[str]:
    """
    Compute stacking labels ('A', 'B', 'C', ...) for atomic layers.

    Implements: MATH_SPEC §9: Gram matrix normal equations in exact Fraction arithmetic,
    followed by floor-modulo reduction s_k^mod = s_k - floor(s_k) in [0, 1).
    Units: Positions and vectors in units of lattice parameter a.
    Exact: Exact rational arithmetic via fractions.Fraction.

    Args:
        layers: Dictionary mapping layer_key -> atom_indices OR list of layer keys / atom indices.
        positions_frac: Atomic positions as Fraction 3-tuples (or NumPy array).
        v1_frac: First 2D primitive in-plane vector as Fraction 3-tuple.
        v2_frac: Second 2D primitive in-plane vector as Fraction 3-tuple.

    Returns:
        List of single-character strings ('A', 'B', 'C', ...) for each layer in ascending d.
    """
    frac_positions = _ensure_fraction_tuples(positions_frac)

    if isinstance(layers, dict):
        sorted_keys = sorted(layers.keys())
        if not sorted_keys:
            return []
        # Per MATH_SPEC §9.1: r0 is any atom in layer 0 (or closest to 0)
        ref_key = 0 if 0 in layers else min(sorted_keys, key=lambda k: abs(k))
        r0 = frac_positions[layers[ref_key][0]]
        rep_atoms = [frac_positions[layers[k][0]] for k in sorted_keys]
    else:
        # layers is a list of atom indices
        if len(layers) == 0:
            return []
        if len(frac_positions) == len(layers):
            rep_atoms = frac_positions
        else:
            rep_atoms = [frac_positions[idx] for idx in layers]
        sorted_keys = list(range(len(rep_atoms)))
        r0 = rep_atoms[0]

    if not rep_atoms:
        return []

    # Gram matrix entries for (v1, v2)
    def dot_f(u: tuple[Fraction, Fraction, Fraction], v: tuple[Fraction, Fraction, Fraction]) -> Fraction:
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

    g11 = dot_f(v1_frac, v1_frac)
    g12 = dot_f(v1_frac, v2_frac)
    g22 = dot_f(v2_frac, v2_frac)
    det = g11 * g22 - g12 * g12

    if det == Fraction(0, 1):
        raise ValueError("In-plane vectors v1 and v2 must be linearly independent.")

    distinct_offsets: list[tuple[Fraction, Fraction]] = []
    label_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    # Compute offset for each atom
    offsets: list[tuple[Fraction, Fraction]] = []
    for rep in rep_atoms:
        delta = (rep[0] - r0[0], rep[1] - r0[1], rep[2] - r0[2])
        b1 = dot_f(v1_frac, delta)
        b2 = dot_f(v2_frac, delta)

        # Solve Gram normal equations by Cramer's rule
        s1 = (b1 * g22 - b2 * g12) / det
        s2 = (g11 * b2 - g12 * b1) / det

        # Exact floor-modulo reduction into [0, 1)
        s1_mod = s1 - Fraction(math.floor(s1))
        s2_mod = s2 - Fraction(math.floor(s2))

        offsets.append((s1_mod, s2_mod))

    # Per MATH_SPEC §9.3: Labels A, B, C... are registered in order of appearance
    # starting from layer 0 upwards (d >= 0)
    order_indices = [i for i, k in enumerate(sorted_keys) if k >= 0]
    neg_indices = [i for i, k in enumerate(sorted_keys) if k < 0]
    registration_order = order_indices + list(reversed(neg_indices))

    for idx in registration_order:
        off = offsets[idx]
        if off not in distinct_offsets:
            distinct_offsets.append(off)

    # Now assign labels
    labels: list[str] = []
    for off in offsets:
        matched_idx = distinct_offsets.index(off)
        labels.append(label_chars[matched_idx % len(label_chars)])

    return labels
