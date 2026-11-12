"""
Unit tests for core/density.py (TP-30 through TP-47).
"""
from __future__ import annotations

import math
import pytest

from lattice_generator.core.density import (
    formula_spacing,
    true_spacing,
    planar_density,
)


def test_tp30_formula_spacing_100():
    """TP-30 formula_spacing(1,0,0) = 1.0 (a=1)."""
    assert math.isclose(formula_spacing(1, 0, 0), 1.0, rel_tol=1e-12)


def test_tp31_formula_spacing_110():
    """TP-31 formula_spacing(1,1,0) = 1/sqrt(2)."""
    assert math.isclose(formula_spacing(1, 1, 0), 1.0 / math.sqrt(2.0), rel_tol=1e-12)


def test_tp32_formula_spacing_111():
    """TP-32 formula_spacing(1,1,1) = 1/sqrt(3)."""
    assert math.isclose(formula_spacing(1, 1, 1), 1.0 / math.sqrt(3.0), rel_tol=1e-12)


def test_tp33_true_spacing_fcc_111():
    """TP-33 true_spacing('fcc', 1,1,1) = 1/sqrt(3) [no extinction: all odd]."""
    assert math.isclose(true_spacing("fcc", 1, 1, 1), 1.0 / math.sqrt(3.0), rel_tol=1e-12)


def test_tp34_true_spacing_fcc_110():
    """TP-34 true_spacing('fcc', 1,1,0) = 1/(2*sqrt(2)) [extinction: mixed parity]."""
    assert math.isclose(true_spacing("fcc", 1, 1, 0), 1.0 / (2.0 * math.sqrt(2.0)), rel_tol=1e-12)


def test_tp35_true_spacing_fcc_200():
    """TP-35 true_spacing('fcc', 2,0,0) = 0.5 [no extinction: all even]."""
    assert math.isclose(true_spacing("fcc", 2, 0, 0), 0.5, rel_tol=1e-12)


def test_tp36_true_spacing_fcc_100():
    """TP-36 true_spacing('fcc', 1,0,0) = 0.5 [extinction: mixed parity]."""
    assert math.isclose(true_spacing("fcc", 1, 0, 0), 0.5, rel_tol=1e-12)


def test_tp37_true_spacing_bcc_110():
    """TP-37 true_spacing('bcc', 1,1,0) = 1/sqrt(2) [h+k+l=2 even, no extinction]."""
    assert math.isclose(true_spacing("bcc", 1, 1, 0), 1.0 / math.sqrt(2.0), rel_tol=1e-12)


def test_tp38_true_spacing_bcc_100():
    """TP-38 true_spacing('bcc', 1,0,0) = 0.5 [h+k+l=1 odd, extinction]."""
    assert math.isclose(true_spacing("bcc", 1, 0, 0), 0.5, rel_tol=1e-12)


def test_tp39_true_spacing_bcc_111():
    """TP-39 true_spacing('bcc', 1,1,1) = 1/(2*sqrt(3)) [h+k+l=3 odd, extinction]."""
    assert math.isclose(true_spacing("bcc", 1, 1, 1), 1.0 / (2.0 * math.sqrt(3.0)), rel_tol=1e-12)


def test_tp40_planar_density_fcc_111():
    """TP-40 planar_density('fcc', 1,1,1) = 4/sqrt(3) ≈ 2.309401... (tol 1e-10)."""
    assert math.isclose(planar_density("fcc", 1, 1, 1), 4.0 / math.sqrt(3.0), rel_tol=1e-10)


def test_tp41_planar_density_fcc_100():
    """TP-41 planar_density('fcc', 1,0,0) = 2.0 (exact)."""
    assert math.isclose(planar_density("fcc", 1, 0, 0), 2.0, rel_tol=1e-10)


def test_tp42_planar_density_fcc_110():
    """TP-42 planar_density('fcc', 1,1,0) = sqrt(2) ≈ 1.414213... (tol 1e-10)."""
    assert math.isclose(planar_density("fcc", 1, 1, 0), math.sqrt(2.0), rel_tol=1e-10)


def test_tp43_planar_density_bcc_110():
    """TP-43 planar_density('bcc', 1,1,0) = sqrt(2) ≈ 1.414213... (tol 1e-10)."""
    assert math.isclose(planar_density("bcc", 1, 1, 0), math.sqrt(2.0), rel_tol=1e-10)


def test_tp44_planar_density_bcc_100():
    """TP-44 planar_density('bcc', 1,0,0) = 1.0 (exact)."""
    assert math.isclose(planar_density("bcc", 1, 0, 0), 1.0, rel_tol=1e-10)


def test_tp45_planar_density_bcc_111():
    """TP-45 planar_density('bcc', 1,1,1) = 1/sqrt(3) ≈ 0.577350... (tol 1e-10)."""
    assert math.isclose(planar_density("bcc", 1, 1, 1), 1.0 / math.sqrt(3.0), rel_tol=1e-10)


def test_tp46_planar_density_antipodal():
    """TP-46 planar_density(-h,-k,-l) = planar_density(h,k,l) for all test cases."""
    for structure in ["fcc", "bcc"]:
        for h, k, l in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0)]:
            pos_dens = planar_density(structure, h, k, l)
            neg_dens = planar_density(structure, -h, -k, -l)
            assert math.isclose(pos_dens, neg_dens, rel_tol=1e-12)


def test_tp47_planar_density_equivalent_family():
    """TP-47 Equivalent planes (100), (010), (001) give identical planar density."""
    for structure in ["fcc", "bcc"]:
        d100 = planar_density(structure, 1, 0, 0)
        d010 = planar_density(structure, 0, 1, 0)
        d001 = planar_density(structure, 0, 0, 1)
        assert math.isclose(d100, d010, rel_tol=1e-12)
        assert math.isclose(d100, d001, rel_tol=1e-12)
