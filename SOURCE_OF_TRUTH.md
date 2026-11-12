# SOURCE OF TRUTH — Lattice Generator

> **This is the single authoritative anchor for the entire project.**
> All other documentation defers to this file. If any other document contradicts
> this one, this document is correct. Update this file at the start of every phase.

---

## 1. Project Scope

**Lattice Generator** is a Python mini-project that:

1. Projects atomic arrangements of **FCC** and **BCC** crystals onto arbitrary {hkl}
   planes.
2. Builds and analyses **bicrystals** (grain-boundary junctions) using the same
   projection engine.
3. Provides a **graphical GUI** as a thin wrapper over the pure-math core.

**Non-negotiable priority:** mathematical correctness. Speed, aesthetics, and
convenience are secondary. If they conflict with correctness, correctness wins.
Test tolerances are never loosened to make a test pass; the implementation is
fixed instead.

---

## 2. Fundamental Conventions

| Symbol | Meaning |
|--------|---------|
| `a` | Lattice parameter (all internal lengths in units of `a`) |
| `a*` | Reciprocal lattice vector `2π/a x̂` (for cubic: `a* = b* = c* = 2π/a`) |
| `hkl` | Miller indices, reduced by gcd, stored as plain Python `int` |
| `n̂` | Unit normal to (hkl) plane |
| `d` | Signed distance from the (hkl) reference plane (atom at origin) |
| `e1, e2` | Orthonormal in-plane basis vectors (right-handed with `n̂`) |

### 2.1 Crystal Structures

**FCC** — conventional cubic cell, 4 atoms per cell:

| Atom | Fractional coordinates |
|------|----------------------|
| 0 | (0, 0, 0) |
| 1 | (½, ½, 0) |
| 2 | (½, 0, ½) |
| 3 | (0, ½, ½) |

**BCC** — conventional cubic cell, 2 atoms per cell:

| Atom | Fractional coordinates |
|------|----------------------|
| 0 | (0, 0, 0) |
| 1 | (½, ½, ½) |

### 2.2 Arithmetic Policy

| Quantity | Required type |
|----------|--------------|
| Plane membership, layer index | Exact (`fractions.Fraction` or doubled-integer coords) |
| Layer assignment, stacking offsets | Exact |
| Sigma, coincidence | Exact integer arithmetic |
| Rendered 2D positions, bond lengths for display | `float` (cast only at render boundary) |
| Rotation matrix entries for CSL presets | Exact rational (`fractions.Fraction`) |
| Rotation matrix entries for arbitrary angles | `float` (numpy), but `R^T R = I` and `det R = +1` asserted |

**Never** use a float tolerance to group atoms into layers where an exact
method (integer/rational arithmetic) exists.

### 2.3 Coordinate System

- **Right-handed** Cartesian frame throughout: `x̂ = [1,0,0]`, `ŷ = [0,1,0]`,
  `ẑ = [0,0,1]`.
- **Conventional cubic cell**: Cartesian axes aligned with lattice vectors
  `a1 = a x̂`, `a2 = a ŷ`, `a3 = a ẑ`. Internally `a = 1`.
- For general (non-cubic) lattices the reciprocal basis is used explicitly;
  the cubic case is validated as a special case.

---

## 3. Module Layout (canonical)

```
lattice_generator/
├── core/            # Pure math — ZERO GUI or matplotlib imports
│   ├── __init__.py
│   ├── lattice.py          # Bravais basis, supercell generation
│   ├── plane.py            # hkl validation, normal, reciprocal lattice
│   ├── projection.py       # Signed distance, orthonormal basis, 2D coords
│   ├── layers.py           # Layer grouping (exact), stacking sequence
│   ├── cell2d.py           # 2D primitive cell, Gauss reduction
│   ├── density.py          # Interplanar spacing, planar density
│   ├── neighbors.py        # Nearest-neighbor distance and CN
│   └── bicrystal.py        # Rotation, CSL/Sigma, bicrystal construction
├── viz/             # Matplotlib/plot code only
│   ├── __init__.py
│   ├── plot2d.py
│   ├── plot_bicrystal.py
│   └── plot_comparison.py
├── export/          # I/O only
│   ├── __init__.py
│   └── exporters.py        # PNG/SVG/CSV/XYZ/JSON
├── gui/             # Thin GUI wrapper — no math
│   ├── __init__.py
│   └── app.py              # PySide6 application
└── tests/
    ├── __init__.py
    ├── test_plane.py
    ├── test_projection.py
    ├── test_layers.py
    ├── test_cell2d.py
    ├── test_density.py
    ├── test_neighbors.py
    ├── test_bicrystal.py
    └── test_integration.py
```

---

## 4. Phase Status Log

| Phase | Description | Status | Approval |
|-------|-------------|--------|----------|
| 0 | Documentation | ✅ Complete | Complete |
| 1 | core/ library + tests | ✅ Complete | Gate Passed (54/54 passed) |
| 2 | Bicrystal engine + tests | ✅ Complete | Gate Passed (16/16 passed) |
| 3 | viz/ and export/ | ✅ Complete | Gate Passed (5/5 passed) |
| 4 | gui/ thin layer | ✅ Complete | Gate Passed |
| 5 | Final audit + README | ✅ Complete | Complete |

### Phase 1 Execution Log
- **Date**: 2026-11-10
- **Tests run**: `pytest tests/ -v --ignore=tests/test_bicrystal.py --ignore=tests/test_integration.py`
- **Results**: 54 passed, 0 failed, 0 skipped.
- **Components**: `core/lattice.py`, `core/plane.py`, `core/projection.py`, `core/layers.py`, `core/cell2d.py`, `core/density.py`, `core/neighbors.py`.
- **Issues resolved**: Fixed boundary truncation in finite supercell layer enumeration by evaluating bulk sequence from reference plane $d=0$ upwards. Verified exact rational consistency $\rho^2 A^2 = 1$ and all property-based tests.

### Phase 2 Execution Log
- **Date**: 2026-11-10
- **Tests run**: `pytest tests/test_bicrystal.py -v`
- **Results**: 16 passed, 0 failed, 0 skipped. Total suite passing: 70/70.
- **Components**: `core/bicrystal.py` (Rodrigues rotations, CSL presets with exact Fraction entries, integer Sigma algorithm with 2-adic removal, bicrystal construction, overlap removal, boundary classification).
- **Issues resolved**: Ensured 0° rotation sets all atoms into single crystal grain with 0 removals. Verified cubic symmetry 90° rotation invariance on FCC lattice.

### Phase 3 Execution Log
- **Date**: 2026-11-10
- **Tests run**: `pytest tests/test_integration.py -v`
- **Results**: 5 passed, 0 failed, 0 skipped. Total suite passing: 75/75.
- **Components**: `viz/plot_2d.py`, `viz/plot_3d.py`, `viz/palettes.py`, `export/writer.py` (CSV, XYZ, JSON, PNG, SVG export).
- **Issues resolved**: Implemented custom JSON encoder for Fraction and NumPy types. Verified 1e-10 round-trip position accuracy and exact layer index preservation across CSV, XYZ, and JSON formats. Validated PNG signature and SVG XML schema.

### Phase 4 Execution Log
- **Date**: 2026-11-10
- **Components**: `gui/app.py`, `gui/__init__.py`, `gui/__main__.py`, `tests/test_gui.py`.
- **Verification**: Complete thin-shell GUI architecture connecting core calculation engines and matplotlib canvas backends. Implemented single-crystal and bicrystal misorientation modes, layer filtering, interactive metrics display, and file export triggers. Headless tests verified under `QT_QPA_PLATFORM=offscreen`.

### Phase 5 Execution Log
- **Date**: 2026-11-10
- **Components**: `README.md`, `SOURCE_OF_TRUTH.md`, code cleanliness audit.
- **Verification**: Completed comprehensive README with mathematical foundations, analytical reference tables, installation instructions, GUI launch guidance, and API usage examples. Audited zero matplotlib/PySide6 imports in `core/`. All non-negotiable mathematical conventions verified.

---

## 5. Key Decisions (with rationale)

| Decision | Rationale |
|----------|-----------|
| GUI framework: **PySide6** | Qt is mature, offline, natively cross-platform. Avoids web-server complexity of Streamlit. See ARCHITECTURE.md §5 for full comparison. |
| Exact arithmetic for layers/sigma | Floating-point grouping would produce wrong layer counts for high-index planes. Exact arithmetic is mandatory per the project brief. |
| `fractions.Fraction` for CSL rotations | Ensures `Sigma·R` is exactly integer with no rounding. |
| Units of `a` internally | Keeps all formulas dimensionless; only the GUI/export layer inserts physical `a`. |
| Gauss reduction for 2D cell | Finds the provably shortest pair of basis vectors; deterministic and exact. |
| Core has zero GUI/matplotlib imports | Testability and separation of concerns. |

---

## 6. Reference Values (ground truth for tests)

### Planar Density (units: `1/a²`)

| Structure | Plane | Value | Exact form |
|-----------|-------|-------|------------|
| FCC | (111) | 2.3094… | 4/√3 |
| FCC | (100) | 2.0 | 2 |
| FCC | (110) | 1.4142… | √2 |
| BCC | (110) | 1.4142… | √2 |
| BCC | (100) | 1.0 | 1 |
| BCC | (111) | 0.5774… | 1/√3 |

### Interplanar Spacing (units: `a`)

Defined as the **actual** separation between consecutive distinct atomic planes
(accounting for extinction):

| Structure | Plane | d_formula | d_actual | Extinction? |
|-----------|-------|-----------|----------|-------------|
| FCC | (111) | a/√3 | a/(2√3) = a/√12 | Yes — mixed parity does not apply, (111) is all-odd → no extinction; actual = a/√3 |
| FCC | (200) | a/2 | a/2 | No |
| FCC | (110) | a/√2 | a/(2√2) | Yes |
| BCC | (110) | a/√2 | a/√2 | No (h+k+l=2, even) |
| BCC | (200) | a/2 | a/2 | No |
| BCC | (111) | a/√3 | a/(2√3) | Yes |

> **Extinction rule (cubic):**
> - FCC: planes with mixed parity (hkl not all-even or all-odd) have zero structure
>   factor; the true inter-*atomic*-plane spacing is half the formula spacing.
> - BCC: planes with h+k+l odd have zero structure factor; true spacing is half.

### Nearest-Neighbour Distances

| Structure | NN distance | CN |
|-----------|------------|-----|
| FCC | a/√2 | 12 |
| BCC | a√3/2 | 8 |

### Stacking Sequences

| Structure | Plane | Sequence |
|-----------|-------|---------|
| FCC | (111) | ABC |
| FCC | (100) | AB |
| FCC | (110) | AB |
| BCC | (110) | AB |
| BCC | (100) | AB |
| BCC | (111) | ABC |

### Bicrystal Sigma Values

| Axis | Angle | Sigma |
|------|-------|-------|
| [100] | 36.87° | 5 |
| [111] | 60° | 3 |

---

## 7. External Dependencies

| Package | Purpose | Required in |
|---------|---------|------------|
| `numpy` | Float linear algebra, rotation matrices | core/, viz/ |
| `fractions` | Exact rational arithmetic (stdlib) | core/ |
| `math` | Integer sqrt, gcd (stdlib) | core/ |
| `matplotlib` | 2D/3D plotting | viz/ |
| `PySide6` | GUI | gui/ |
| `pytest` | Test runner | tests/ |
| `hypothesis` | Property-based tests | tests/ |

---

## 8. Glossary

| Term | Definition |
|------|-----------|
| **CSL** | Coincidence Site Lattice: the set of lattice sites shared by both grains after rotation |
| **Sigma (Σ)** | Ratio of CSL unit-cell volume to crystal unit-cell volume; smallest positive integer such that Σ·R has integer entries (after removing common factors; odd part for cubic) |
| **Extinction** | A reflection (hkl) for which the structure factor is zero; atoms on intermediate planes destructively interfere |
| **Gauss reduction** | Algorithm that finds the shortest two independent vectors spanning a 2D lattice |
| **Tilt boundary** | GB where the rotation axis lies in the boundary plane (axis · n̂ = 0) |
| **Twist boundary** | GB where the rotation axis is parallel to the boundary normal (axis ∥ n̂) |
| **Planar density** | Number of atoms per unit area in one atomic plane (units: 1/a²) |
| **Layer** | Set of all atoms sharing the same exact signed distance from the reference (hkl) plane |
