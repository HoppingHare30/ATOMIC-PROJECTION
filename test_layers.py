"""
Unit and property-based tests for core/layers.py (TP-13 through TP-23, TP-PB-4, TP-PB-7).
"""
from __future__ import annotations

from fractions import Fraction
import math
import numpy as np
import pytest
from hypothesis import given, settings
import hypothesis.strategies as st

from lattice_generator.core.lattice import get_basis, generate_supercell_exact
from lattice_generator.core.layers import (
    assign_layers,
    filter_by_n_layers,
    filter_by_distance,
    stacking_sequence,
)
from lattice_generator.core.cell2d import find_inplane_vectors


def canonical_relabel(seq: list[str]) -> list[str]:
    """Relabel sequence to A, B, C... in order of first appearance."""
    mapping: dict[str, str] = {}
    result: list[str] = []
    chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for item in seq:
        if item not in mapping:
            mapping[item] = chars[len(mapping)]
        result.append(mapping[item])
    return result


def test_tp13_assign_layers_fcc_111():
    """TP-13 assign_layers(FCC, 1,1,1, 3x3x3 supercell): layer keys are exact integers with uniform step."""
    atoms = generate_supercell_exact("fcc", 3)
    layer_dict = assign_layers(atoms, 1, 1, 1)

    sorted_keys = sorted(layer_dict.keys())
    assert len(sorted_keys) > 0
    assert all(isinstance(k, int) for k in sorted_keys)

    # Check that keys have regular integer step Delta_d_int = 2 (no missing layers)
    diffs = [sorted_keys[i + 1] - sorted_keys[i] for i in range(len(sorted_keys) - 1)]
    assert all(d == 2 for d in diffs)


def test_tp14_fcc_111_2x2x2_distinct_layers():
    """TP-14 FCC (111) 2x2x2 block has 6 distinct layers with step a/sqrt(3)."""
    # 2x2x2 conventional unit cell block: coordinates in [0, 2]^3 with d < 2*sqrt(3)
    basis = get_basis("fcc")
    atoms_2x2x2 = []
    for n1 in [0, 1, 2]:
        for n2 in [0, 1, 2]:
            for n3 in [0, 1, 2]:
                for tau in basis:
                    x = Fraction(n1) + tau[0]
                    y = Fraction(n2) + tau[1]
                    z = Fraction(n3) + tau[2]
                    if x <= 2 and y <= 2 and z <= 2:
                        atoms_2x2x2.append((x, y, z))

    # Deduplicate atoms
    atoms_2x2x2 = list(set(atoms_2x2x2))
    layer_dict = assign_layers(atoms_2x2x2, 1, 1, 1)

    # 6 distinct layers from d=0 to d < 2*sqrt(3) across the 2x2x2 cell
    norm = math.sqrt(3.0)
    pos_keys = [k for k in sorted(layer_dict.keys()) if 0 <= k < 12]
    assert len(pos_keys) == 6

    # Spacing between consecutive layers: delta_d = delta_key / (2 * sqrt(3))
    diffs = [pos_keys[i + 1] - pos_keys[i] for i in range(len(pos_keys) - 1)]
    assert all(d == 2 for d in diffs)
    d_step = diffs[0] / (2.0 * norm)
    assert math.isclose(d_step, 1.0 / norm, rel_tol=1e-12)


def test_tp15_filter_by_n_layers():
    """TP-15 filter_by_n_layers(n_layers=3): at most 3 distinct layers returned per side."""
    atoms = generate_supercell_exact("fcc", 4)
    layer_dict = assign_layers(atoms, 1, 1, 1)
    filtered = filter_by_n_layers(layer_dict, atoms, 1, 1, 1, n_layers=3)

    keys = sorted(filtered.keys())
    neg_or_zero = [k for k in keys if k <= 0]
    pos = [k for k in keys if k > 0]
    assert len(neg_or_zero) == 3
    assert len(pos) == 3


def test_tp16_filter_by_distance():
    """TP-16 filter_by_distance range [0, 0.5]: all returned atoms satisfy 0 <= d <= 0.5a."""
    atoms = generate_supercell_exact("fcc", 3)
    layer_dict = assign_layers(atoms, 1, 1, 1)
    filtered = filter_by_distance(layer_dict, atoms, 1, 1, 1, d_min=0.0, d_max=0.5)

    norm = math.sqrt(3.0)
    for key, indices in filtered.items():
        d = key / (2.0 * norm)
        assert -1e-12 <= d <= 0.5 + 1e-12
        for idx in indices:
            pos = atoms[idx]
            d_atom = float(pos[0] + pos[1] + pos[2]) / norm
            assert -1e-12 <= d_atom <= 0.5 + 1e-12


def _get_stacking_labels(structure: str, h: int, k: int, l: int, N: int = 3) -> list[str]:
    """Helper to extract repeating stacking sequence for a given crystal plane."""
    atoms = generate_supercell_exact(structure, N)
    layer_dict = assign_layers(atoms, h, k, l)
    # Take contiguous bulk layers starting from key 0 upwards
    sorted_pos_keys = [k for k in sorted(layer_dict.keys()) if k >= 0]
    step = sorted_pos_keys[1] - sorted_pos_keys[0] if len(sorted_pos_keys) > 1 else 1
    contiguous = [sorted_pos_keys[0]]
    for k_val in sorted_pos_keys[1:]:
        if k_val - contiguous[-1] == step:
            contiguous.append(k_val)
        else:
            break
    bulk_layers = {k: layer_dict[k] for k in contiguous}
    v1, v2 = find_inplane_vectors(structure, h, k, l)
    return stacking_sequence(bulk_layers, atoms, v1, v2)


def test_tp17_stacking_fcc_111():
    """TP-17 stacking_sequence(FCC, 1,1,1) has period-3 'ABC' repeat."""
    labels = _get_stacking_labels("fcc", 1, 1, 1, N=4)
    # Check that sequence contains A, B, C repeating
    assert len(set(labels)) == 3
    # First 3 labels are A, B, C
    assert labels[:3] == ["A", "B", "C"]
    # Periodicity check
    for i in range(len(labels) - 3):
        assert labels[i] == labels[i + 3]


def test_tp18_stacking_fcc_100():
    """TP-18 stacking_sequence(FCC, 1,0,0) has period-2 'AB' repeat."""
    labels = _get_stacking_labels("fcc", 1, 0, 0, N=3)
    assert len(set(labels)) == 2
    assert labels[:2] == ["A", "B"]
    for i in range(len(labels) - 2):
        assert labels[i] == labels[i + 2]


def test_tp19_stacking_fcc_110():
    """TP-19 stacking_sequence(FCC, 1,1,0) has period-2 'AB' repeat."""
    labels = _get_stacking_labels("fcc", 1, 1, 0, N=3)
    assert len(set(labels)) == 2
    assert labels[:2] == ["A", "B"]
    for i in range(len(labels) - 2):
        assert labels[i] == labels[i + 2]


def test_tp20_stacking_bcc_110():
    """TP-20 stacking_sequence(BCC, 1,1,0) has period-2 'AB' repeat."""
    labels = _get_stacking_labels("bcc", 1, 1, 0, N=3)
    assert len(set(labels)) == 2
    assert labels[:2] == ["A", "B"]
    for i in range(len(labels) - 2):
        assert labels[i] == labels[i + 2]


def test_tp21_stacking_bcc_100():
    """TP-21 stacking_sequence(BCC, 1,0,0) has period-2 'AB' repeat."""
    labels = _get_stacking_labels("bcc", 1, 0, 0, N=3)
    assert len(set(labels)) == 2
    assert labels[:2] == ["A", "B"]
    for i in range(len(labels) - 2):
        assert labels[i] == labels[i + 2]


def test_tp22_stacking_bcc_111():
    """TP-22 stacking_sequence(BCC, 1,1,1) has period-3 'ABC' repeat."""
    labels = _get_stacking_labels("bcc", 1, 1, 1, N=4)
    assert len(set(labels)) == 3
    assert labels[:3] == ["A", "B", "C"]
    for i in range(len(labels) - 3):
        assert labels[i] == labels[i + 3]


def test_tp23_antipodal_distances():
    """TP-23 (-h,-k,-l) gives same set of distinct layer distances (magnitudes) as (h,k,l)."""
    for structure in ["fcc", "bcc"]:
        atoms = generate_supercell_exact(structure, 2)
        for h, k, l in [(1, 1, 1), (1, 1, 0), (2, 1, 0)]:
            ld_pos = assign_layers(atoms, h, k, l)
            ld_neg = assign_layers(atoms, -h, -k, -l)

            norm = math.sqrt(h * h + k * k + l * l)
            d_pos = sorted(set(abs(k / (2.0 * norm)) for k in ld_pos.keys()))
            d_neg = sorted(set(abs(k / (2.0 * norm)) for k in ld_neg.keys()))
            assert np.allclose(d_pos, d_neg, atol=1e-12)


@st.composite
def hkl_strategy(draw):
    h = draw(st.integers(min_value=-5, max_value=5))
    k = draw(st.integers(min_value=-5, max_value=5))
    l = draw(st.integers(min_value=-5, max_value=5))
    if h == 0 and k == 0 and l == 0:
        h = 1
    return (h, k, l)


@settings(max_examples=100)
@given(hkl=hkl_strategy())
def test_tp_pb4_layer_exact_property(hkl):
    """TP-PB-4 Layer assignment produces exact integer keys for all atoms."""
    h, k, l = hkl
    for structure in ["fcc", "bcc"]:
        atoms = generate_supercell_exact(structure, 1)
        layer_dict = assign_layers(atoms, h, k, l)
        assert all(isinstance(k_val, int) for k_val in layer_dict.keys())


def test_tp_pb7_antipodal_stacking_invariant():
    """TP-PB-7 Antipodal invariant: period length equal and canonical_relabel(reversed(seq_pos)) == seq_neg."""
    test_cases = [
        ("fcc", 1, 1, 1),
        ("fcc", 1, 1, 0),
        ("fcc", 1, 0, 0),
        ("bcc", 1, 1, 0),
        ("bcc", 1, 0, 0),
        ("bcc", 1, 1, 1),
    ]
    for structure, h, k, l in test_cases:
        seq_pos = _get_stacking_labels(structure, h, k, l, N=4)
        seq_neg = _get_stacking_labels(structure, -h, -k, -l, N=4)

        # Unique period length
        p_pos = len(set(seq_pos))
        p_neg = len(set(seq_neg))
        assert p_pos == p_neg

        # Canonical relabelling of the reversal of one period
        period_pos = seq_pos[:p_pos]
        period_neg = seq_neg[:p_neg]
        rev_relabelled = canonical_relabel(list(reversed(period_pos)))
        assert rev_relabelled == period_neg
