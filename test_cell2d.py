"""
Unit and property-based tests for core/cell2d.py (TP-24 through TP-29, TP-PB-5).
"""
from __future__ import annotations

from fractions import Fraction
import math
import numpy as np
import pytest
from hypothesis import given, settings
import hypothesis.strategies as st

from lattice_generator.core.cell2d import (
    find_inplane_vectors,
    gauss_reduce,
    cell_area,
    cell_area_sq,
    _cross_frac,
    _dot_frac,
)
from lattice_generator.core.density import (
    true_spacing,
    true_spacing_sq,
    planar_density,
    planar_density_sq,
)


def test_tp24_parity_and_perpendicularity():
    """TP-24 In-plane vectors satisfy perpendicularity and structure parity rules."""
    test_cases = [
        ("fcc", 1, 1, 1),
        ("fcc", 1, 1, 0),
        ("fcc", 1, 0, 0),
        ("bcc", 1, 1, 0),
        ("bcc", 1, 0, 0),
        ("bcc", 1, 1, 1),
        ("sc", 1, 0, 0),
    ]
    for structure, h, k, l in test_cases:
        v1, v2 = find_inplane_vectors(structure, h, k, l)
        for v in (v1, v2):
            w1 = int(2 * v[0])
            w2 = int(2 * v[1])
            w3 = int(2 * v[2])

            # Integer perpendicularity
            assert h * w1 + k * w2 + l * w3 == 0

            # Parity check
            if structure == "sc":
                assert w1 % 2 == 0 and w2 % 2 == 0 and w3 % 2 == 0
            elif structure == "fcc":
                assert (w1 + w2 + w3) % 2 == 0
            elif structure == "bcc":
                assert w1 % 2 == w2 % 2 == w3 % 2


def test_tp25_linear_independence():
    """TP-25 In-plane vectors are linearly independent (cross product != 0 in Fraction)."""
    for structure in ["fcc", "bcc", "sc"]:
        for h, k, l in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]:
            v1, v2 = find_inplane_vectors(structure, h, k, l)
            cross = _cross_frac(v1, v2)
            assert cross != (Fraction(0), Fraction(0), Fraction(0))


def test_tp26_gauss_reduce_minkowski():
    """TP-26 gauss_reduce returns |v1|² <= |v2|² and (v1·v2)/|v1|² in (-0.5, 0.5]."""
    test_vectors = [
        (
            (Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),
            (Fraction(3, 1), Fraction(2, 1), Fraction(0, 1)),
        ),
        (
            (Fraction(1, 2), Fraction(1, 2), Fraction(0, 1)),
            (Fraction(1, 1), Fraction(0, 1), Fraction(1, 2)),
        ),
        (
            (Fraction(2, 1), Fraction(1, 1), Fraction(0, 1)),
            (Fraction(1, 1), Fraction(2, 1), Fraction(0, 1)),
        ),
    ]
    for v1_in, v2_in in test_vectors:
        v1, v2 = gauss_reduce(v1_in, v2_in)
        len1 = _dot_frac(v1, v1)
        len2 = _dot_frac(v2, v2)
        assert len1 <= len2

        ratio = _dot_frac(v1, v2) / len1
        assert Fraction(-1, 2) < ratio <= Fraction(1, 2)


def test_tp27_cell_area_fcc_111():
    """TP-27 cell_area(FCC, 1,1,1) = sqrt(3)/4 ≈ 0.4330127 a² (tol 1e-10)."""
    v1, v2 = find_inplane_vectors("fcc", 1, 1, 1)
    area = cell_area(v1, v2)
    expected = math.sqrt(3.0) / 4.0
    assert math.isclose(area, expected, rel_tol=1e-10)


def test_tp28_tp_pb5_rho_area_consistency():
    """TP-28 / TP-PB-5 Verify rho^2 * A_2D^2 == Fraction(1, 1) and math.isclose(rho * A_2D, 1.0)."""
    test_planes = [
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 1, 0), (1, 0, 1), (0, 1, 1),
        (1, 1, 1), (-1, 1, 1), (1, -1, 1),
        (2, 1, 0), (1, 2, 0), (2, 1, 1),
        (3, 1, 0), (3, 1, 1),
    ]
    for structure in ["fcc", "bcc", "sc"]:
        n_conv = 4 if structure == "fcc" else (2 if structure == "bcc" else 1)
        for h, k, l in test_planes:
            v1, v2 = find_inplane_vectors(structure, h, k, l)
            a_sq = cell_area_sq(v1, v2)
            d_sq = true_spacing_sq(structure, h, k, l)
            rho_sq = Fraction(n_conv * n_conv, 1) * d_sq

            # Non-negotiable exact rational consistency check (MATH_SPEC §7.4)
            assert rho_sq * a_sq == Fraction(1, 1), f"Failed for {structure} ({h},{k},{l}): {rho_sq} * {a_sq} != 1"

            # Float consistency sanity check
            a_float = cell_area(v1, v2)
            rho_float = planar_density(structure, h, k, l)
            assert math.isclose(rho_float * a_float, 1.0, rel_tol=1e-12)


def test_tp29_equivalent_planes_area():
    """TP-29 Equivalent planes under cubic symmetry give identical cell area."""
    for structure in ["fcc", "bcc"]:
        planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0)]
        areas = []
        for h, k, l in planes:
            v1, v2 = find_inplane_vectors(structure, h, k, l)
            areas.append(cell_area_sq(v1, v2))
        assert all(a == areas[0] for a in areas)

        planes110 = [(1, 1, 0), (1, 0, 1), (0, 1, 1), (-1, 1, 0)]
        areas110 = []
        for h, k, l in planes110:
            v1, v2 = find_inplane_vectors(structure, h, k, l)
            areas110.append(cell_area_sq(v1, v2))
        assert all(a == areas110[0] for a in areas110)
