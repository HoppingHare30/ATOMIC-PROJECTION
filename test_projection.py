"""
Unit and property-based tests for core/projection.py (TP-08 through TP-12, TP-PB-2, TP-PB-3).
"""
from __future__ import annotations

from fractions import Fraction
import numpy as np
import pytest
from hypothesis import given, settings
import hypothesis.strategies as st

from lattice_generator.core.lattice import get_basis, generate_supercell
from lattice_generator.core.plane import unit_normal
from lattice_generator.core.projection import orthonormal_basis, project_2d, signed_distance_exact


@st.composite
def hkl_strategy(draw):
    h = draw(st.integers(min_value=-5, max_value=5))
    k = draw(st.integers(min_value=-5, max_value=5))
    l = draw(st.integers(min_value=-5, max_value=5))
    if h == 0 and k == 0 and l == 0:
        h = 1
    return (h, k, l)


def test_tp08_orthonormal_basis_111():
    """TP-08 orthonormal_basis(1,1,1): |e1|=|e2|=1, e1·e2=0, e1·n=0, e2·n=0 to 1e-12."""
    e1, e2, n = orthonormal_basis(1, 1, 1)
    assert abs(np.linalg.norm(e1) - 1.0) < 1e-12
    assert abs(np.linalg.norm(e2) - 1.0) < 1e-12
    assert abs(np.dot(e1, e2)) < 1e-12
    assert abs(np.dot(e1, n)) < 1e-12
    assert abs(np.dot(e2, n)) < 1e-12


def test_tp09_right_handedness():
    """TP-09 e1 x e2 = n_hat for (1,1,1), (1,0,0), (1,1,0)."""
    for h, k, l in [(1, 1, 1), (1, 0, 0), (1, 1, 0)]:
        e1, e2, n = orthonormal_basis(h, k, l)
        cross = np.cross(e1, e2)
        assert np.allclose(cross, n, atol=1e-12)


@settings(max_examples=100)
@given(hkl=hkl_strategy())
def test_tp10_tp_pb2_orthonormal_property(hkl):
    """TP-10 / TP-PB-2 For 100 random (h,k,l): all orthonormality conditions hold to 1e-12."""
    h, k, l = hkl
    e1, e2, n = orthonormal_basis(h, k, l)
    assert abs(np.linalg.norm(e1) - 1.0) < 1e-12
    assert abs(np.linalg.norm(e2) - 1.0) < 1e-12
    assert abs(np.linalg.norm(n) - 1.0) < 1e-12
    assert abs(np.dot(e1, e2)) < 1e-12
    assert abs(np.dot(e1, n)) < 1e-12
    assert abs(np.dot(e2, n)) < 1e-12
    assert np.allclose(np.cross(e1, e2), n, atol=1e-12)


def test_tp11_tp_pb3_reconstruction():
    """TP-11 / TP-PB-3 3D reconstruction proj_x*e1 + proj_y*e2 + d*n_hat matches original to 1e-12."""
    for structure in ["fcc", "bcc"]:
        positions = generate_supercell(structure, 2)
        for h, k, l in [(1, 1, 1), (1, 1, 0), (1, 0, 0), (2, 1, 0)]:
            e1, e2, n = orthonormal_basis(h, k, l)
            proj = project_2d(positions, e1, e2)
            d = np.dot(positions, n)

            # Reconstruct 3D position
            recon = (
                proj[:, 0:1] * e1.reshape(1, 3)
                + proj[:, 1:2] * e2.reshape(1, 3)
                + d[:, np.newaxis] * n.reshape(1, 3)
            )
            assert np.allclose(recon, positions, atol=1e-12)

            # Invariant 2: n_hat . r_perp_3D = 0 to 1e-12
            r_perp_3d = proj[:, 0:1] * e1.reshape(1, 3) + proj[:, 1:2] * e2.reshape(1, 3)
            dot_residuals = np.dot(r_perp_3d, n)
            assert np.all(np.abs(dot_residuals) < 1e-12)


def test_tp12_signed_distance_exact_integer():
    """TP-12 signed_distance_exact: 2*(hx+ky+lz) is an exact integer for all basis atoms."""
    for structure in ["fcc", "bcc"]:
        basis = get_basis(structure)
        for h, k, l in [(1, 1, 1), (1, 1, 0), (1, 0, 0), (3, 2, 1), (-1, 2, 4)]:
            for atom in basis:
                d_prime = signed_distance_exact(atom, h, k, l)
                key_frac = Fraction(2, 1) * d_prime
                assert key_frac.denominator == 1
                assert isinstance(int(key_frac), int)
