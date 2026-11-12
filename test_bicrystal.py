"""
Tests for Bicrystal Construction, Rodrigues Rotations, and CSL Presets.

Implements unit tests TP-53 through TP-67 and property-based test TP-PB-6
per TEST_PLAN.md §3 and §6.
"""
from __future__ import annotations

from fractions import Fraction
import math
from hypothesis import given, settings
import hypothesis.strategies as st
import numpy as np
import pytest

from lattice_generator.core.bicrystal import (
    build_bicrystal,
    classify_boundary,
    coincidence_cell_area,
    coincident_site_fraction,
    compute_sigma,
    extract_rotation_axis_angle,
    get_csl_preset,
    rodrigues_rotation,
)
from lattice_generator.core.cell2d import cell_area, find_inplane_vectors
from lattice_generator.core.neighbors import nearest_neighbour


def test_tp53_rotation_matrix_0_deg():
    """TP-53 rotation_matrix([0,0,1], 0°) -> identity matrix to 1e-12."""
    R = rodrigues_rotation([0.0, 0.0, 1.0], 0.0)
    assert np.allclose(R, np.eye(3), atol=1e-12)


def test_tp54_rotation_matrix_90_deg():
    """TP-54 rotation_matrix([1,0,0], 90°): R^T R = I, det R = +1 [tol 1e-12]."""
    R = rodrigues_rotation([1.0, 0.0, 0.0], math.pi / 2.0)
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert math.isclose(np.linalg.det(R), 1.0, abs_tol=1e-12)


def test_tp55_csl_rotation_sigma5():
    """TP-55 csl_rotation('sigma5_001') and 'sigma5_100': R^T R = I, det R = +1 [exact Fraction]."""
    # sigma5_001: axis is [001]
    R_001 = get_csl_preset("sigma5_001")
    assert R_001.dtype == object

    # Exact R^T R = I
    prod_001 = R_001.T @ R_001
    for i in range(3):
        for j in range(3):
            expected = Fraction(1, 1) if i == j else Fraction(0, 1)
            assert prod_001[i, j] == expected

    # Exact det R = 1
    # det(R) = R00*(R11*R22 - R12*R21) - R01*(R10*R22 - R12*R20) + R02*(R10*R21 - R11*R20)
    det_001 = (
        R_001[0, 0] * (R_001[1, 1] * R_001[2, 2] - R_001[1, 2] * R_001[2, 1])
        - R_001[0, 1] * (R_001[1, 0] * R_001[2, 2] - R_001[1, 2] * R_001[2, 0])
        + R_001[0, 2] * (R_001[1, 0] * R_001[2, 1] - R_001[1, 1] * R_001[2, 0])
    )
    assert det_001 == Fraction(1, 1)

    # Fixed axis is [001]: R * [0, 0, 1]^T == [0, 0, 1]^T
    v_z = np.array([Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)], dtype=object)
    assert np.array_equal(R_001 @ v_z, v_z)

    # sigma5_100: axis is [100]
    R_100 = get_csl_preset("sigma5_100")
    prod_100 = R_100.T @ R_100
    for i in range(3):
        for j in range(3):
            expected = Fraction(1, 1) if i == j else Fraction(0, 1)
            assert prod_100[i, j] == expected

    v_x = np.array([Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)], dtype=object)
    assert np.array_equal(R_100 @ v_x, v_x)


def test_tp56_csl_rotation_sigma3():
    """TP-56 csl_rotation('sigma3_111'): R^T R = I, det R = +1, Tr = 2 -> 60° [exact Fraction]."""
    R = get_csl_preset("sigma3_111")
    assert R.dtype == object

    # Exact orthogonality
    prod = R.T @ R
    for i in range(3):
        for j in range(3):
            expected = Fraction(1, 1) if i == j else Fraction(0, 1)
            assert prod[i, j] == expected

    # Exact determinant
    det = (
        R[0, 0] * (R[1, 1] * R[2, 2] - R[1, 2] * R[2, 1])
        - R[0, 1] * (R[1, 0] * R[2, 2] - R[1, 2] * R[2, 0])
        + R[0, 2] * (R[1, 0] * R[2, 1] - R[1, 1] * R[2, 0])
    )
    assert det == Fraction(1, 1)

    # Trace is 2 -> Tr = 1 + 2*cos(theta) -> cos(theta) = 1/2 -> theta = 60°
    tr = R[0, 0] + R[1, 1] + R[2, 2]
    assert tr == Fraction(2, 1)

    # Fixed axis is [111]
    v_111 = np.array([Fraction(1, 1), Fraction(1, 1), Fraction(1, 1)], dtype=object)
    assert np.array_equal(R @ v_111, v_111)


def test_tp57_compute_sigma_sigma5():
    """TP-57 compute_sigma('sigma5_001') and 'sigma5_100' both equal 5 [exact integer]."""
    R_001 = get_csl_preset("sigma5_001")
    assert compute_sigma(R_001) == 5

    R_100 = get_csl_preset("sigma5_100")
    assert compute_sigma(R_100) == 5


def test_tp58_compute_sigma_sigma3():
    """TP-58 compute_sigma('sigma3_111') = 3 [exact integer]."""
    R_111 = get_csl_preset("sigma3_111")
    assert compute_sigma(R_111) == 3


def test_tp59_build_bicrystal_0_deg():
    """TP-59 build_bicrystal(0° rotation): all atoms in one grain, removed_count=0."""
    R_0 = rodrigues_rotation([0.0, 0.0, 1.0], 0.0)
    res = build_bicrystal("fcc", 1, 1, 1, R_0, N=2)

    assert res.removed_count == 0
    assert len(res.grain2_positions) == 0
    assert len(res.grain1_positions) == len(res.positions)
    assert np.all(res.grain_ids == 1)


def test_tp60_build_bicrystal_90_deg_cubic_symmetry():
    """TP-60 build_bicrystal(90° about [100]): result identical to single crystal."""
    R_90 = rodrigues_rotation([1.0, 0.0, 0.0], math.pi / 2.0)
    res = build_bicrystal("fcc", 1, 0, 0, R_90, N=3)

    assert res.removed_count == 0

    # In FCC, every atom r must have 2r as integers with sum-even
    w = 2.0 * res.positions
    w_round = np.round(w)
    assert np.allclose(w, w_round, atol=1e-10)
    parities = np.sum(w_round.astype(int), axis=1) % 2
    assert np.all(parities == 0)

    # Nearest neighbour distance across interface matches FCC d_NN
    d_nn, _ = nearest_neighbour("fcc")
    min_interface_dist = float("inf")
    for p2 in res.grain2_positions:
        d = np.min(np.linalg.norm(res.grain1_positions - p2, axis=1))
        if d < min_interface_dist:
            min_interface_dist = d
    assert math.isclose(min_interface_dist, d_nn, rel_tol=1e-10)


def test_tp61_grain_distances():
    """TP-61 Grain 1 atoms all have d <= 0; Grain 2 atoms all have d > 0."""
    R = get_csl_preset("sigma5_001")
    res = build_bicrystal("fcc", 0, 0, 1, R, N=2)

    # Normal is [0, 0, 1]
    assert len(res.grain1_positions) > 0
    assert len(res.grain2_positions) > 0

    # Grain 1: z <= 0
    assert np.all(res.grain1_positions[:, 2] <= 1e-12)

    # Grain 2: z > 0
    assert np.all(res.grain2_positions[:, 2] > 1e-12)


def test_tp62_removal_threshold():
    """TP-62 Removal threshold = 0.5 * d_NN; no kept atom has d < 0.5*d_NN to Grain 1."""
    R = rodrigues_rotation([0.0, 0.0, 1.0], math.radians(20.0))
    res = build_bicrystal("fcc", 0, 0, 1, R, N=2, removal_alpha=0.5)

    d_nn, _ = nearest_neighbour("fcc")
    thresh = 0.5 * d_nn

    # All kept Grain 2 atoms must have distance >= threshold to every Grain 1 atom
    if len(res.grain2_positions) > 0:
        for p2 in res.grain2_positions:
            min_d = np.min(np.linalg.norm(res.grain1_positions - p2, axis=1))
            assert min_d >= thresh - 1e-10


def test_tp63_boundary_classification():
    """TP-63 Boundary classification: axis perp n -> 'tilt'; axis || n -> 'twist'."""
    normal = np.array([0.0, 0.0, 1.0])

    axis_tilt = np.array([0.0, 1.0, 0.0])
    assert classify_boundary(axis_tilt, normal) == "tilt"

    axis_twist = np.array([0.0, 0.0, 1.0])
    assert classify_boundary(axis_twist, normal) == "twist"

    axis_mixed = np.array([1.0, 1.0, 1.0])
    assert classify_boundary(axis_mixed, normal) == "mixed"


def test_tp64_sigma_invariance():
    """TP-64 Sigma invariant: equivalent axis/angle representations give same Sigma."""
    # Sigma 5 can be represented as 36.87° about [001] or [100]
    R_001 = get_csl_preset("sigma5_001")
    R_100 = get_csl_preset("sigma5_100")
    assert compute_sigma(R_001) == compute_sigma(R_100) == 5

    # Also via Rodrigues float representation
    theta_rad = math.atan2(3.0, 4.0)
    R_rodrigues = rodrigues_rotation([0.0, 0.0, 1.0], theta_rad)
    assert compute_sigma(R_rodrigues) == 5


def test_tp65_coincidence_cell_area():
    """TP-65 Coincidence-cell area in plan view = Sigma * A_2D(reference grain)."""
    # For FCC (111): A_2D = sqrt(3)/4
    v1, v2 = find_inplane_vectors("fcc", 1, 1, 1)
    a_2d = cell_area(v1, v2)
    assert math.isclose(a_2d, math.sqrt(3.0) / 4.0, rel_tol=1e-12)

    sigma = 3
    a_csl = coincidence_cell_area(sigma, a_2d)
    assert math.isclose(a_csl, 3.0 * (math.sqrt(3.0) / 4.0), rel_tol=1e-12)


def test_tp66_fraction_of_coincident_sites():
    """TP-66 Fraction of coincident sites in plan view = 1 / Sigma."""
    assert math.isclose(coincident_site_fraction(3), 1.0 / 3.0, rel_tol=1e-12)
    assert math.isclose(coincident_site_fraction(5), 1.0 / 5.0, rel_tol=1e-12)
    assert math.isclose(coincident_site_fraction(7), 1.0 / 7.0, rel_tol=1e-12)


def test_tp67_misorientation_angle():
    """TP-67 Misorientation angle extracted from R matches input angle to 0.01°."""
    test_angles = [5.0, 15.2, 36.87, 45.0, 60.0, 70.53, 90.0]
    for angle_deg in test_angles:
        R = rodrigues_rotation([1.0, 2.0, 3.0], math.radians(angle_deg))
        axis, theta_rad = extract_rotation_axis_angle(R)
        extracted_deg = math.degrees(theta_rad)
        assert abs(extracted_deg - angle_deg) < 0.01


@settings(max_examples=100)
@given(
    axis=st.tuples(
        st.floats(min_value=-10.0, max_value=10.0),
        st.floats(min_value=-10.0, max_value=10.0),
        st.floats(min_value=-10.0, max_value=10.0),
    ).filter(lambda v: abs(v[0]) + abs(v[1]) + abs(v[2]) > 0.1),
    theta=st.floats(min_value=-2.0 * math.pi, max_value=2.0 * math.pi),
)
def test_tp_pb6_rodrigues_property(axis, theta):
    """TP-PB-6 Random axis + angle: R^T R = I and det R = +1 [tol 1e-12]."""
    R = rodrigues_rotation(axis, theta)
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert math.isclose(np.linalg.det(R), 1.0, abs_tol=1e-12)
