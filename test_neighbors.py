"""
Unit tests for core/neighbors.py (TP-48 through TP-52).
"""
from __future__ import annotations

from fractions import Fraction
import math
import pytest

from lattice_generator.core.neighbors import nearest_neighbour


def test_tp48_nearest_neighbour_fcc_distance():
    """TP-48 nearest_neighbour('fcc')[0] = 1/sqrt(2) ≈ 0.7071067... (tol 1e-10)."""
    d_nn, _ = nearest_neighbour("fcc")
    expected = 1.0 / math.sqrt(2.0)
    assert math.isclose(d_nn, expected, rel_tol=1e-10)


def test_tp49_nearest_neighbour_fcc_cn():
    """TP-49 nearest_neighbour('fcc')[1] = 12 (exact integer)."""
    _, cn = nearest_neighbour("fcc")
    assert cn == 12
    assert isinstance(cn, int)


def test_tp50_nearest_neighbour_bcc_distance():
    """TP-50 nearest_neighbour('bcc')[0] = sqrt(3)/2 ≈ 0.8660254... (tol 1e-10)."""
    d_nn, _ = nearest_neighbour("bcc")
    expected = math.sqrt(3.0) / 2.0
    assert math.isclose(d_nn, expected, rel_tol=1e-10)


def test_tp51_nearest_neighbour_bcc_cn():
    """TP-51 nearest_neighbour('bcc')[1] = 8 (exact integer)."""
    _, cn = nearest_neighbour("bcc")
    assert cn == 8
    assert isinstance(cn, int)


def test_tp52_nearest_neighbour_adapts_to_custom_geometry():
    """TP-52 d_NN computed from geometry, not hardcoded (verified on simple cubic basis)."""
    # Simple cubic basis has NN distance = 1.0 and CN = 6 (along +/-x, +/-y, +/-z)
    sc_basis = [(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1))]
    d_sc, cn_sc = nearest_neighbour(sc_basis)
    assert math.isclose(d_sc, 1.0, rel_tol=1e-10)
    assert cn_sc == 6
