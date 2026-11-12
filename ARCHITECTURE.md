# ARCHITECTURE — Lattice Generator

> This document describes the module layout, data flow, interfaces, and
> key architectural decisions. It defers to SOURCE_OF_TRUTH.md for conventions
> and to MATH_SPEC.md for mathematical definitions.

---

## §1 Guiding Principles

1. **Correctness first.** The `core/` package contains zero GUI or plotting
   imports; it can be tested in a headless environment with no display.
2. **Single projection engine.** Features 1–3, 7, 8, 10, and the entire
   bicrystal builder all call the same `core/projection.py` routines. No
   duplicate geometry code.
3. **Exact arithmetic at the boundary.** Layer grouping, stacking-offset
   comparison, Sigma computation, and CSL classification use
   `fractions.Fraction` or plain integers. The `float` cast happens only when
   passing data to `viz/` or `export/`.
4. **GUI is a thin shell.** `gui/` calls `core/` and passes results to `viz/`
   and `export/`. No formulae, no geometry, no business logic in `gui/`.
5. **Tests exercise `core/` directly.** Tests import from `core/` only; they
   do not test GUI or plot code (those are integration-tested separately).

---

## §2 Package Layout

```
lattice_generator/
│
├── core/                   # Pure math — no GUI, no matplotlib
│   ├── __init__.py
│   ├── lattice.py          # Atom position generation
│   ├── plane.py            # hkl validation, GCD reduction, normal vector
│   ├── projection.py       # Signed distance, orthonormal basis, 2D coords
│   ├── layers.py           # Exact layer grouping, stacking sequence
│   ├── cell2d.py           # In-plane lattice vectors, Gauss reduction
│   ├── density.py          # d_hkl (formula + measured), planar density
│   ├── neighbors.py        # NN distance and coordination number
│   └── bicrystal.py        # Rotation, CSL, Sigma, bicrystal construction
│
├── viz/                    # Matplotlib code only
│   ├── __init__.py
│   ├── plot2d.py           # 2D scatter, layer colours, 2D cell overlay
│   ├── plot_bicrystal.py   # Plan view and side view of bicrystal
│   └── plot_comparison.py  # Side-by-side multi-plane comparison
│
├── export/                 # I/O only, no math
│   ├── __init__.py
│   └── exporters.py        # PNG/SVG via matplotlib; CSV, XYZ, JSON writers
│
├── gui/                    # PySide6 application — no math
│   ├── __init__.py
│   └── app.py
│
└── tests/
    ├── __init__.py
    ├── test_plane.py        # §3.1–3.3 of MATH_SPEC
    ├── test_projection.py   # §6 of MATH_SPEC
    ├── test_layers.py       # §4, §9 of MATH_SPEC
    ├── test_cell2d.py       # §7 of MATH_SPEC
    ├── test_density.py      # §5, §8 of MATH_SPEC
    ├── test_neighbors.py    # §10 of MATH_SPEC
    ├── test_bicrystal.py    # §11–§13 of MATH_SPEC
    └── test_integration.py  # Cross-module consistency checks
```

---

## §3 Module Interfaces

### §3.1 `core/lattice.py`

```python
def get_basis(structure: str) -> list[tuple[Fraction, Fraction, Fraction]]:
    """
    Return the fractional-coordinate basis atoms for the given structure.
    structure: 'fcc' or 'bcc'.
    Returns list of (u, v, w) with Fraction entries.
    MATH_SPEC §1.2.
    """

def generate_supercell(
    structure: str,
    n1: int, n2: int, n3: int
) -> np.ndarray:  # shape (N, 3), dtype float64, units of a
    """
    Generate atom positions for an n1×n2×n3 supercell.
    Positions in units of a. MATH_SPEC §1.3.
    """
```

### §3.2 `core/plane.py`

```python
def validate_and_reduce(h: int, k: int, l: int) -> tuple[int, int, int]:
    """
    Validate (hkl) != (000); reduce by gcd(|h|,|k|,|l|).
    Returns reduced (h, k, l). MATH_SPEC §3.3.
    """

def unit_normal(h: int, k: int, l: int) -> np.ndarray:  # shape (3,)
    """
    Unit normal n = (h,k,l)/sqrt(h²+k²+l²). MATH_SPEC §3.2.
    """
```

### §3.3 `core/projection.py`

```python
def signed_distance_exact(
    pos: tuple[Fraction, Fraction, Fraction],
    h: int, k: int, l: int
) -> Fraction:
    """
    Exact signed distance (without normalisation) = h*x + k*y + l*z.
    Returns a Fraction proportional to d*sqrt(h²+k²+l²). MATH_SPEC §4.
    """

def orthonormal_basis(
    h: int, k: int, l: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns (e1, e2, n_hat) — orthonormal right-handed frame. MATH_SPEC §6.1.
    """

def project_2d(
    positions: np.ndarray,  # (N, 3) in units of a
    e1: np.ndarray,
    e2: np.ndarray
) -> np.ndarray:  # (N, 2)
    """
    Project 3D positions onto the plane via e1, e2. MATH_SPEC §6.2.
    """
```

### §3.4 `core/layers.py`

```python
def assign_layers(
    positions_frac: list[tuple[Fraction, Fraction, Fraction]],
    h: int, k: int, l: int
) -> dict[int, list[int]]:
    """
    Group atom indices by exact layer key = 2*(h*x + k*y + l*z).
    Returns dict mapping integer layer key → list of atom indices.
    MATH_SPEC §4 (exact arithmetic policy).
    """

def stacking_sequence(
    layers: dict[int, list[int]],
    positions_2d: np.ndarray,
    v1_2d: np.ndarray,
    v2_2d: np.ndarray
) -> list[str]:
    """
    Compute stacking label (A/B/C/…) per layer. MATH_SPEC §9.
    """
```

### §3.5 `core/cell2d.py`

```python
def find_inplane_vectors(h: int, k: int, l: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Find the two shortest independent in-plane lattice vectors by enumeration
    + Gauss reduction. MATH_SPEC §7.
    """

def gauss_reduce(v1: np.ndarray, v2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Gauss/Lagrange reduction of a 2D lattice basis. MATH_SPEC §7.2.
    """

def cell_area(v1: np.ndarray, v2: np.ndarray) -> float:
    """
    Area of the 2D cell |v1 × v2|. MATH_SPEC §7.3.
    """
```

### §3.6 `core/density.py`

```python
def formula_spacing(h: int, k: int, l: int) -> float:
    """
    d_hkl = a/sqrt(h²+k²+l²) with a=1. MATH_SPEC §5.1.
    """

def true_spacing(
    structure: str,
    h: int, k: int, l: int
) -> float:
    """
    True interplanar spacing from atomic layer positions. MATH_SPEC §5.2–5.3.
    """

def planar_density(
    structure: str,
    h: int, k: int, l: int,
    n_supercell: int = 10
) -> float:
    """
    Atoms per unit area in one (hkl) plane. MATH_SPEC §8.
    """
```

### §3.7 `core/neighbors.py`

```python
def nearest_neighbour(structure: str) -> tuple[float, int]:
    """
    Returns (d_NN, CN) computed from geometry. MATH_SPEC §10.2.
    Not hardcoded: computed from supercell distances.
    """
```

### §3.8 `core/bicrystal.py`

```python
def rotation_matrix(axis: np.ndarray, angle_deg: float) -> np.ndarray:
    """
    Rodrigues rotation matrix. MATH_SPEC §11.
    Float implementation; asserts R^T R = I and det R = +1.
    """

def csl_rotation(preset: str) -> np.ndarray:
    """
    Return exact rational CSL rotation matrix as np.ndarray of Fraction.
    preset: 'sigma5_100' or 'sigma3_111'. MATH_SPEC §11.1.
    """

def compute_sigma(R: np.ndarray) -> int:
    """
    Compute Sigma from rotation matrix. MATH_SPEC §12.2.
    """

def build_bicrystal(
    structure: str,
    h: int, k: int, l: int,
    R: np.ndarray,
    n_layers: int,
    removal_alpha: float = 0.5
) -> dict:
    """
    Construct the bicrystal slab. Returns dict with keys:
      'grain1_positions', 'grain2_positions', 'removed_count',
      'boundary_type', 'sigma', 'misorientation_deg'.
    MATH_SPEC §13.
    """
```

---

## §4 Data Flow

```
User input (hkl, structure, params)
        │
        ▼
  core/plane.py          ← validates hkl, computes n̂
        │
        ▼
  core/lattice.py        ← generates supercell positions
        │
        ├──► core/projection.py    ← (e1, e2, n̂), 2D positions, signed distances
        │           │
        │           ├──► core/layers.py     ← exact layer groups, stacking
        │           ├──► core/cell2d.py     ← 2D primitive cell
        │           └──► core/density.py    ← d_hkl, planar density
        │
        ├──► core/neighbors.py     ← NN distance, CN
        │
        └──► core/bicrystal.py     ← rotation, Sigma, bicrystal slab
                    │
                    └── (reuses projection.py for plan/side views)
        │
        ▼
  viz/ ─────────────────────────────────────► matplotlib figures
  export/ ───────────────────────────────────► PNG/SVG/CSV/XYZ/JSON
        │
        ▼
  gui/ (PySide6) ──────────► widgets, callbacks, display panels
```

---

## §5 GUI Framework: PySide6

### §5.1 Comparison

| Criterion | PySide6 | Streamlit | PyQt5 | tkinter |
|-----------|---------|-----------|-------|---------|
| Offline (no browser/server) | ✅ | ❌ (needs server) | ✅ | ✅ |
| Native look & feel | ✅ | ❌ | ✅ | ⚠️ |
| matplotlib embedding | ✅ (FigureCanvas) | ✅ | ✅ | ✅ |
| LGPL license | ✅ | Apache 2 | GPL | Python |
| Modern widget set | ✅ | ✅ | ✅ | ❌ |
| Complex layouts | ✅ | ❌ | ✅ | ⚠️ |
| Signal/slot architecture | ✅ | ❌ | ✅ | ❌ |

**Decision: PySide6.** It is offline, natively cross-platform, has a mature
layout system for the complex multi-panel interface required, and embeds
matplotlib figures natively via `FigureCanvasQTAgg`. Streamlit's server
model and limited layout control make it unsuitable for a tool requiring
real-time parameter sweeping.

### §5.2 GUI Layout

```
┌─ Lattice Generator ───────────────────────────────────────────────────┐
│ ┌─ Controls ───────────┐  ┌─ Viewer ─────────────────────────────────┐ │
│ │ Structure: [FCC▼]    │  │                                           │ │
│ │ h: [1] k: [1] l: [1]│  │          [matplotlib canvas]              │ │
│ │ Layers: [5]          │  │                                           │ │
│ │ ─────────────────── │  └───────────────────────────────────────────┘ │
│ │ [x] Layer colours    │  ┌─ Metrics ────────────────────────────────┐ │
│ │ [ ] NN bonds         │  │ d_hkl:        0.333 a                    │ │
│ │ [x] Unit cell        │  │ Planar dens:  2.309 /a²                  │ │
│ │ ─────────────────── │  │ Stacking:     ABC                        │ │
│ │ ── Bicrystal ──────  │  │ NN distance:  0.707 a                    │ │
│ │ Axis: [100]          │  └───────────────────────────────────────────┘ │
│ │ Angle: [36.87]       │                                               │
│ │ Preset: [Sigma5▼]   │                                               │
│ │ View: [Plan▼]        │                                               │
│ │ ─────────────────── │                                               │
│ │ [Export PNG] [CSV]   │                                               │
│ └──────────────────────┘                                               │
└───────────────────────────────────────────────────────────────────────┘
```

---

## §6 Dependency Policy

| Import from... | Allowed in `core/`? | Allowed in `viz/`? | Allowed in `gui/`? |
|----------------|--------------------|--------------------|-------------------|
| `numpy` | ✅ | ✅ | ✅ (for display only) |
| `fractions` | ✅ | ❌ (use float) | ❌ |
| `math` | ✅ | ✅ | ✅ |
| `matplotlib` | ❌ | ✅ | ❌ (use `FigureCanvas`) |
| `PySide6` | ❌ | ❌ | ✅ |
| Any `core` module | ✅ | ✅ | ✅ |
| Any `viz` module | ❌ | ✅ | ✅ |
| Any `gui` module | ❌ | ❌ | ✅ |

---

## §7 Testing Strategy (summary; full plan in TEST_PLAN.md)

- `pytest` for all unit and integration tests.
- `hypothesis` for property-based tests (random hkl).
- Tests import from `core/` only.
- Zero float-tolerance loosening: tolerances fixed to `1e-12` (machine epsilon
  for double arithmetic) and never changed to make a test pass.
- Full test suite runs at the end of every phase; results reported before
  the next phase begins.

---

## §8 File and Format Conventions

| Artefact | Format | Location |
|----------|--------|---------|
| 2D plot (raster) | PNG (300 dpi) | export/ |
| 2D plot (vector) | SVG | export/ |
| Atom data | CSV (index, x, y, z, layer, proj_x, proj_y) | export/ |
| 3D coordinates | XYZ (standard format) | export/ |
| Computed metrics | JSON (d, density, stacking, sigma, …) | export/ |

---

## §9 Build and Environment

```
# requirements.txt (minimum)
numpy>=1.26
matplotlib>=3.8
PySide6>=6.6
pytest>=8.0
hypothesis>=6.100
```

No package manager beyond pip is required. The project is a standard Python
package installed with `pip install -e .` from the repo root.
