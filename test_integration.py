"""
Integration and Round-Trip Tests (Phase 3).

Implements TP-68 through TP-72 per TEST_PLAN.md §4.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import pytest

from lattice_generator.core.lattice import generate_supercell
from lattice_generator.core.layers import assign_layers
from lattice_generator.core.projection import orthonormal_basis, project_2d
from lattice_generator.core.plane import unit_normal
from lattice_generator.export.writer import (
    export_csv,
    export_image,
    export_json,
    export_xyz,
)
from lattice_generator.viz.plot_2d import plot_plan_view, plot_side_view


@pytest.fixture
def sample_crystal_data():
    """Fixture providing atomic coordinates, layer keys, and 2D projections."""
    atoms = generate_supercell("fcc", 2)
    h, k, l = 1, 1, 1
    layer_dict = assign_layers(atoms, h, k, l)

    # Flatten layer keys into array aligned with atoms
    layer_keys = np.zeros(len(atoms), dtype=int)
    for key, indices in layer_dict.items():
        for idx in indices:
            layer_keys[idx] = key

    # Projections
    e1, e2, n_hat = orthonormal_basis(h, k, l)
    proj_2d = project_2d(atoms, e1, e2)
    u = proj_2d[:, 0]
    v = proj_2d[:, 1]
    d = np.dot(atoms, n_hat)

    return {
        "atoms": atoms,
        "layer_keys": layer_keys,
        "u": u,
        "v": v,
        "d": d,
        "hkl": (h, k, l),
    }


def test_tp68_export_csv_round_trip(tmp_path: Path, sample_crystal_data):
    """TP-68 Export CSV -> reimport -> positions match to 1e-10, layer indices exact."""
    csv_file = tmp_path / "test_crystal.csv"
    atoms = sample_crystal_data["atoms"]
    layer_keys = sample_crystal_data["layer_keys"]

    export_csv(csv_file, atoms, layer_keys)
    assert csv_file.exists()

    reimported_pos = []
    reimported_keys = []
    with open(csv_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            reimported_pos.append([float(row["x"]), float(row["y"]), float(row["z"])])
            reimported_keys.append(int(row["layer_index"]))

    reimported_pos = np.array(reimported_pos)
    assert len(reimported_pos) == len(atoms)
    assert np.allclose(reimported_pos, atoms, atol=1e-10)
    assert reimported_keys == list(layer_keys)


def test_tp69_export_xyz_round_trip(tmp_path: Path, sample_crystal_data):
    """TP-69 Export XYZ -> reimport -> 3D positions match to 1e-10."""
    xyz_file = tmp_path / "test_crystal.xyz"
    atoms = sample_crystal_data["atoms"]

    export_xyz(xyz_file, atoms, species="Cu", comment="FCC (111) test")
    assert xyz_file.exists()

    reimported_pos = []
    with open(xyz_file, "r", encoding="utf-8") as f:
        lines = f.readlines()
        count = int(lines[0].strip())
        assert count == len(atoms)
        comment = lines[1].strip()
        assert comment == "FCC (111) test"
        for line in lines[2:]:
            parts = line.split()
            if len(parts) == 4:
                species, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
                assert species == "Cu"
                reimported_pos.append([x, y, z])

    reimported_pos = np.array(reimported_pos)
    assert np.allclose(reimported_pos, atoms, atol=1e-10)


def test_tp70_export_json_round_trip(tmp_path: Path, sample_crystal_data):
    """TP-70 Export JSON -> reimport -> all numeric values match to 1e-10."""
    json_file = tmp_path / "test_data.json"
    data = {
        "structure": "fcc",
        "hkl": [1, 1, 1],
        "planar_density": 2.309401076758503,
        "interplanar_spacing": 0.5773502691896258,
        "positions": sample_crystal_data["atoms"][:10].tolist(),
        "layer_keys": sample_crystal_data["layer_keys"][:10].tolist(),
    }

    export_json(json_file, data)
    assert json_file.exists()

    with open(json_file, "r", encoding="utf-8") as f:
        reimported = json.load(f)

    assert reimported["structure"] == "fcc"
    assert reimported["hkl"] == [1, 1, 1]
    assert math.isclose(reimported["planar_density"], data["planar_density"], abs_tol=1e-10)
    assert math.isclose(reimported["interplanar_spacing"], data["interplanar_spacing"], abs_tol=1e-10)
    assert np.allclose(reimported["positions"], data["positions"], atol=1e-10)
    assert reimported["layer_keys"] == data["layer_keys"]


def test_tp71_export_png(tmp_path: Path, sample_crystal_data):
    """TP-71 Export PNG -> file exists and size > 0 bytes; valid image file."""
    png_file = tmp_path / "plot.png"
    u = sample_crystal_data["u"]
    v = sample_crystal_data["v"]
    layer_keys = sample_crystal_data["layer_keys"]

    fig = plot_plan_view(u, v, layer_keys, title="Test Plan View")
    export_image(fig, png_file, dpi=100)

    assert png_file.exists()
    assert png_file.stat().st_size > 1000  # Non-trivial image file

    # Verify PNG signature bytes: \x89PNG\r\n\x1a\n
    with open(png_file, "rb") as f:
        header = f.read(8)
        assert header == b"\x89PNG\r\n\x1a\n"


def test_tp72_export_svg(tmp_path: Path, sample_crystal_data):
    """TP-72 Export SVG -> valid XML (xml.etree.ElementTree parses without error)."""
    svg_file = tmp_path / "plot.svg"
    u = sample_crystal_data["u"]
    d = sample_crystal_data["d"]

    fig = plot_side_view(u, d, title="Test Side View")
    export_image(fig, svg_file)

    assert svg_file.exists()
    assert svg_file.stat().st_size > 500

    # Parse XML structure
    tree = ET.parse(str(svg_file))
    root = tree.getroot()
    # SVG root tag ends with 'svg'
    assert root.tag.endswith("svg")
