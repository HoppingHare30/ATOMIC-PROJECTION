# Lattice Generator: Crystal Plane & Bicrystal Interface Analyzer

A rigorous computational materials science library and desktop application for generating, analyzing, and visualizing crystallographic planes and bicrystal grain boundaries in cubic crystal systems (FCC, BCC, and SC).

---

## 1. Key Capabilities

- **Single Crystal Planar Analysis**:
  - Exact Miller index reduction and validation $\gcd(|h|, |k|, |l|)$.
  - Analytical interplanar spacings: formula spacing $d_{hkl} = a/\sqrt{h^2+k^2+l^2}$ and true crystallographic spacing $d^\text{true}_{hkl}$ with Bravais extinction rules.
  - Exact planar densities $\rho = \frac{N_\text{conv} \Delta d_\text{int}}{a^2 (h^2+k^2+l^2)}$ with certified consistency $\rho^2 A_{2D}^2 = 1$.
  - In-plane primitive unit cell basis vectors $(v_1, v_2)$ via integer null-space search and Minkowski/Gauss 2D lattice reduction.
  - Layer grouping using doubled-integer exact keys $2(hx+ky+lz) \in \mathbb{Z}$ without floating-point comparisons.
  - Periodic stacking label assignment (e.g., ABC for FCC (111), AB for BCC (110), ABCDEF for BCC (111)) referenced to $d=0$ with antipodal sequence invariance.
  - Coordination number and nearest-neighbour distance ($d_\text{NN} = 1/\sqrt{2}$ for FCC, $\sqrt{3}/2$ for BCC).

- **Bicrystal Grain Boundary Engine**:
  - Proper 3D rotations via Rodrigues' rotation formula with certified $R^T R = I$ and $\det R = +1$.
  - Standard CSL (Coincidence Site Lattice) presets stored as exact `fractions.Fraction` matrices:
    - $\Sigma 3$ $[111]$ ($60.00^\circ$)
    - $\Sigma 5$ $[001]$ ($36.87^\circ$)
    - $\Sigma 5$ $[100]$ ($36.87^\circ$)
    - $\Sigma 7$ $[111]$ ($38.21^\circ$)
    - $\Sigma 9$ $[110]$ ($38.94^\circ$)
    - $\Sigma 11$ $[110]$ ($50.48^\circ$)
  - Integer $\Sigma$ algorithm with exact 2-adic removal.
  - Centered supercell slab construction ($d \le 0$ for Grain 1, $d > 0$ for Grain 2).
  - Overlap atom removal with threshold derived dynamically from $0.5 \times d_\text{NN}$.
  - Grain boundary classification: Tilt ($\mathbf{k} \perp \hat{\mathbf{n}}$), Twist ($\mathbf{k} \parallel \hat{\mathbf{n}}$), or Mixed.
  - Coincidence-cell area $A_\text{CSL} = \Sigma A_{2D}$ and coincident site fraction $1/\Sigma$.

- **Visualization & Export**:
  - 2D Plan View with layer/stacking color and marker encoding, primitive cell parallelograms, and CSL supercells.
  - Side View (tilt/edge-on view) showing atomic layers and the interface boundary at $d=0$.
  - 3D Interactive structure rendering with crystallographic boundary plane highlights.
  - Publication-ready export to CSV, standard XYZ molecular coordinates, JSON metadata, PNG, and SVG vector graphics.

---

## 2. Mathematical Reference Values ($a = 1$)

| Structure | Plane $(hkl)$ | $d_{hkl}$ (Formula) | $d^\text{true}_{hkl}$ (True) | Planar Density $\rho$ ($a^{-2}$) | Stacking Sequence |
|-----------|---------------|---------------------|------------------------------|----------------------------------|-------------------|
| FCC | (1 1 1) | $1/\sqrt{3} \approx 0.577350$ | $1/\sqrt{3} \approx 0.577350$ | $4/\sqrt{3} \approx 2.309401$ | ABC |
| FCC | (1 0 0) | $1.000000$ | $1/2 = 0.500000$ | $2.000000$ | AB |
| FCC | (1 1 0) | $1/\sqrt{2} \approx 0.707107$ | $1/(2\sqrt{2}) \approx 0.353553$ | $2\sqrt{2} \approx 2.828427$ | AB |
| BCC | (1 1 0) | $1/\sqrt{2} \approx 0.707107$ | $1/\sqrt{2} \approx 0.707107$ | $\sqrt{2} \approx 1.414214$ | AB |
| BCC | (1 0 0) | $1.000000$ | $1/2 = 0.500000$ | $1.000000$ | AB |
| BCC | (1 1 1) | $1/\sqrt{3} \approx 0.577350$ | $1/(2\sqrt{3}) \approx 0.288675$ | $2\sqrt{3} \approx 3.464102$ | ABCDEF |

---

## 3. Installation & Requirements

- Python $\ge 3.9$ (tested on Python 3.9 and 3.11+)
- NumPy $\ge 1.26$
- Matplotlib $\ge 3.8$
- PySide6 $\ge 6.6$
- Pytest $\ge 8.0$ and Hypothesis $\ge 6.100$

Install package in editable mode:
```bash
pip install -e .
```

---

## 4. Running the Application

### Graphical User Interface (PySide6)
Launch the GUI application:
```bash
python -m lattice_generator.gui
```

### Programmatic Python API
```python
import numpy as np
from lattice_generator.core.lattice import generate_supercell
from lattice_generator.core.density import planar_density, true_interplanar_spacing
from lattice_generator.core.bicrystal import build_bicrystal, get_csl_preset

# Single crystal planar density
rho = planar_density("fcc", 1, 1, 1)
d_true = true_interplanar_spacing("fcc", 1, 1, 1)
print(f"FCC (111): d = {d_true:.4f} a, rho = {rho:.4f} a^-2")

# Bicrystal Sigma 5 interface
R_sigma5 = get_csl_preset("sigma5_001")
res = build_bicrystal("fcc", 0, 0, 1, R_sigma5, N=3)
print(f"Sigma: {res.sigma}, Boundary: {res.boundary_type}, Removed atoms: {res.removed_count}")
```

---

## 5. Running the Test Suite

The comprehensive automated test suite validates exact analytical reference values, extinction rules, and mathematical invariants:

```bash
pytest tests/ -v
```

Test coverage includes:
- `tests/test_plane.py`: Miller reduction and unit normal antipodal invariants.
- `tests/test_projection.py`: Orthonormal basis right-handedness and reconstruction.
- `tests/test_layers.py`: Exact doubled-integer layer keys and bulk stacking sequences.
- `tests/test_cell2d.py`: Minkowski/Gauss reduced 2D basis and exact consistency $\rho^2 A^2 = 1$.
- `tests/test_density.py`: Interplanar spacings, planar densities, and extinction conditions.
- `tests/test_neighbors.py`: Coordination numbers and nearest-neighbour distances.
- `tests/test_bicrystal.py`: Rodrigues rotations, exact CSL Fraction matrices, Sigma algorithm, overlap removal, and boundary classification.
- `tests/test_integration.py`: High-precision round-trip serialization (CSV, XYZ, JSON, PNG, SVG).
- `tests/test_gui.py`: Headless PySide6 window execution and state updates.
