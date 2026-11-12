"""
Unit and property-based tests for core/plane.py (TP-01 through TP-07, TP-PB-1).
"""
from __future__ import annotations

import math
import numpy as np
import pytest
from hypothesis import given, settings
import hypothesis.strategies as st

from lattice_generator.core.plane import unit_normal, validate_and_reduce


# Hypothesis strategy for valid non-zero (h, k, l) triples
@st.composite
def hkl_strategy(draw):
    h = draw(st.integers(min_value=-5, max_value=5))
    k = draw(st.integers(min_value=-5, max_value=5))
    l = draw(st.integers(min_value=-5, max_value=5))
    if h == 0 and k == 0 and l == 0:
        h = 1
    return (h, k, l)


def test_tp01_validate_and_reduce_zero():
    """TP-01 validate_and_reduce(0, 0, 0) raises ValueError."""
    with pytest.raises(ValueError, match="cannot all be zero"):
        validate_and_reduce(0, 0, 0)


def test_tp02_validate_and_reduce_gcd():
    """TP-02 validate_and_reduce(2, 4, 2) -> (1, 2, 1)."""
    assert validate_and_reduce(2, 4, 2) == (1, 2, 1)


def test_tp03_validate_and_reduce_negatives():
    """TP-03 validate_and_reduce(-1, -1, -1) -> (-1, -1, -1)."""
    assert validate_and_reduce(-1, -1, -1) == (-1, -1, -1)


def test_tp04_unit_normal_111():
    """TP-04 unit_normal(1, 1, 1) -> (1,1,1)/sqrt(3) with norm=1 to 1e-12."""
    n = unit_normal(1, 1, 1)
    expected = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    assert np.allclose(n, expected, atol=1e-12)
    assert abs(np.linalg.norm(n) - 1.0) < 1e-12


def test_tp05_unit_normal_parallel():
    """TP-05 unit_normal(h, k, l) . (h, k, l) > 0."""
    for h, k, l in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (-1, 2, -3)]:
        n = unit_normal(h, k, l)
        dot = np.dot(n, np.array([h, k, l], dtype=np.float64))
        assert dot > 0.0


def test_tp06_unit_normal_antipodal():
    """TP-06 unit_normal(-h, -k, -l) = -unit_normal(h, k, l)."""
    for h, k, l in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 3)]:
        n_pos = unit_normal(h, k, l)
        n_neg = unit_normal(-h, -k, -l)
        assert np.allclose(n_neg, -n_pos, atol=1e-12)


@settings(max_examples=100)
@given(hkl=hkl_strategy())
def test_tp07_tp_pb1_unit_normal_property(hkl):
    """TP-07 / TP-PB-1 |unit_normal(h, k, l)| = 1.0 to 1e-12 for 100 random (h, k, l)."""
    h, k, l = hkl
    n = unit_normal(h, k, l)
    assert abs(np.linalg.norm(n) - 1.0) < 1e-12
    dot = np.dot(n, np.array([h, k, l], dtype=np.float64))
    assert dot > 0.0
