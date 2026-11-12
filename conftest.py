"""
Pytest configuration and environment setup.
"""
from __future__ import annotations

import os
import sys

# Ensure project root and package root are in sys.path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
pkg_root = os.path.join(project_root, "lattice_generator")

if project_root not in sys.path:
    sys.path.insert(0, project_root)
if pkg_root not in sys.path:
    sys.path.insert(0, pkg_root)

# Alias `core` to `lattice_generator.core` so both import styles work
import lattice_generator.core as core_pkg
sys.modules.setdefault("core", core_pkg)
import lattice_generator.core.lattice as lattice_mod
sys.modules.setdefault("core.lattice", lattice_mod)
import lattice_generator.core.plane as plane_mod
sys.modules.setdefault("core.plane", plane_mod)
import lattice_generator.core.projection as proj_mod
sys.modules.setdefault("core.projection", proj_mod)
import lattice_generator.core.layers as layers_mod
sys.modules.setdefault("core.layers", layers_mod)
import lattice_generator.core.cell2d as cell2d_mod
sys.modules.setdefault("core.cell2d", cell2d_mod)
import lattice_generator.core.density as density_mod
sys.modules.setdefault("core.density", density_mod)
import lattice_generator.core.neighbors as neighbors_mod
sys.modules.setdefault("core.neighbors", neighbors_mod)
