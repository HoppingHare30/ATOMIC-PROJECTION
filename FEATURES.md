# FEATURES — Lattice Generator

> Each feature has a canonical identifier (F-01 … F-13), a plain-English
> description, and **acceptance criteria** that map directly to tests in
> TEST_PLAN.md. A feature is "done" only when all its acceptance criteria
> pass. Source-of-truth for conventions: SOURCE_OF_TRUTH.md.
> Source-of-truth for mathematics: MATH_SPEC.md.

---

## F-01 — Plane Definition

**Description:** Accept (hkl) integers from the user; validate, reduce by gcd,
and compute the plane normal.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F01-AC1 | Input (0, 0, 0) raises a `ValueError` with a descriptive message. |
| F01-AC2 | Input (2, 4, 2) is reduced to (1, 2, 1). |
| F01-AC3 | Input (-1, -1, -1) is stored as (-1, -1, -1) (negatives are valid). |
| F01-AC4 | `unit_normal(h, k, l)` has magnitude 1.0 to within 1e-12. |
| F01-AC5 | `unit_normal(h, k, l)` is parallel to `(h, k, l)`. |
| F01-AC6 | `unit_normal(h, k, l)` and `unit_normal(-h, -k, -l)` are antiparallel (differ by sign). |

**Implemented in:** `core/plane.py`  
**Tested in:** `tests/test_plane.py`

---

## F-02 — Orthogonal Projection

**Description:** Project all atoms in a slab onto the (hkl) plane; construct the
orthonormal in-plane basis (e1, e2, n̂) and 2D coordinate array.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F02-AC1 | `|e1| = 1`, `|e2| = 1` to within 1e-12. |
| F02-AC2 | `e1 · e2 = 0` to within 1e-12. |
| F02-AC3 | `e1 · n̂ = 0`, `e2 · n̂ = 0` to within 1e-12. |
| F02-AC4 | `e1 × e2 = n̂` (right-handedness) to within 1e-12. |
| F02-AC5 | For every projected atom, its 3D residual (r − proj) is parallel to n̂. |
| F02-AC6 | The signed distance of a projected point from the reference plane is 0. |
| F02-AC7 | Holds for at least 100 random valid (hkl) inputs (property-based). |

**Implemented in:** `core/projection.py`  
**Tested in:** `tests/test_projection.py`

---

## F-03 — Layer Colouring

**Description:** Group atoms by exact signed distance; assign an integer layer
index; colour layers distinctly in plots.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F03-AC1 | Layer grouping uses exact arithmetic (no float tolerance). |
| F03-AC2 | All atoms within one layer share the same exact signed-distance key. |
| F03-AC3 | FCC (111) supercell yields the correct number of distinct layers for a given slab. |
| F03-AC4 | Layer indices are contiguous integers (no gaps). |
| F03-AC5 | (hkl) and (-h,-k,-l) produce the same set of distinct layers (just reversed sign). |

**Implemented in:** `core/layers.py`  
**Tested in:** `tests/test_layers.py`

---

## F-04 — 2D Plotting

**Description:** Render a scatter plot of projected atom positions, coloured by
layer index, with equal aspect ratio.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F04-AC1 | Plot has equal aspect ratio (`ax.set_aspect('equal')`). |
| F04-AC2 | Each distinct layer is rendered in a visually distinct colour. |
| F04-AC3 | Plot is generated without error for FCC/BCC and any valid (hkl). |
| F04-AC4 | Axes are labelled with physical units (multiples of `a`). |

**Implemented in:** `viz/plot2d.py`  
**Tested in:** `tests/test_integration.py` (smoke test: no exception raised)

---

## F-05 — Slab Thickness Control

**Description:** Allow the user to filter atoms by number of layers (from each
side of the reference plane) or by a distance range along n̂.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F05-AC1 | Requesting N layers returns exactly N distinct layers on each side. |
| F05-AC2 | Requesting a distance range `[d_min, d_max]` returns all atoms with `d_min ≤ d ≤ d_max`. |
| F05-AC3 | An empty slab (no atoms in range) raises a clear `ValueError`. |
| F05-AC4 | Filter is applied before projection (not after). |

**Implemented in:** `core/layers.py`  
**Tested in:** `tests/test_layers.py`

---

## F-06 — Interplanar Spacing and Planar Density

**Description:** Compute d_hkl both from the formula and from atomic layer
separations; assert agreement. Compute planar density.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F06-AC1 | `formula_spacing(h,k,l)` matches `a/sqrt(h²+k²+l²)` to 1e-12. |
| F06-AC2 | `true_spacing(structure, h,k,l)` matches the formula value for non-extinct planes. |
| F06-AC3 | `true_spacing` is half the formula value for FCC mixed-parity planes. |
| F06-AC4 | `true_spacing` is half the formula value for BCC planes with h+k+l odd. |
| F06-AC5 | Planar density matches reference values in SOURCE_OF_TRUTH.md §6 to 1e-10. |
| F06-AC6 | `planar_density × A_2D` is a positive integer for all tested planes. |
| F06-AC7 | (hkl) and (-h,-k,-l) give identical planar density and spacing. |

**Implemented in:** `core/density.py`  
**Tested in:** `tests/test_density.py`

---

## F-07 — 2D Unit Cell Overlay

**Description:** Find the two shortest independent in-plane lattice vectors,
draw the 2D primitive cell on the plot, and report area and angle.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F07-AC1 | The two vectors returned satisfy `h*u + k*v + l*w = 0` (lie in the plane). |
| F07-AC2 | No shorter independent pair exists (Gauss-reduced). |
| F07-AC3 | The cell area matches `A_2D` from §7.3 of MATH_SPEC to 1e-10. |
| F07-AC4 | The vectors are linearly independent (non-zero cross product). |
| F07-AC5 | Cell is drawn on the plot with correct parallelogram shape. |
| F07-AC6 | The angle between the two vectors is reported in degrees. |

**Implemented in:** `core/cell2d.py` + `viz/plot2d.py`  
**Tested in:** `tests/test_cell2d.py`

---

## F-08 — Stacking Sequence Identification

**Description:** For each layer, compute the in-plane offset modulo the 2D
lattice and assign A/B/C labels.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F08-AC1 | FCC (111) gives stacking sequence ABC. |
| F08-AC2 | FCC (100) gives stacking sequence AB. |
| F08-AC3 | FCC (110) gives stacking sequence AB. |
| F08-AC4 | BCC (110) gives stacking sequence AB. |
| F08-AC5 | BCC (100) gives stacking sequence AB. |
| F08-AC6 | BCC (111) gives stacking sequence ABC. |
| F08-AC7 | In-plane offsets are computed with exact arithmetic. |
| F08-AC8 | The sequence repeats with the correct period. |

**Implemented in:** `core/layers.py`  
**Tested in:** `tests/test_layers.py`

---

## F-09 — Multi-Plane Comparison

**Description:** Tabulate and plot side-by-side: d_hkl, planar density, 2D cell
area, stacking sequence, nearest in-plane distance — for a user-specified list
of planes.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F09-AC1 | Table contains all 5 columns (d, density, area, stacking, nearest in-plane). |
| F09-AC2 | Values match per-feature results to 1e-10. |
| F09-AC3 | Side-by-side plots have equal scale (same axis limits if requested). |
| F09-AC4 | At least 3 planes can be compared simultaneously. |

**Implemented in:** `viz/plot_comparison.py` + `core/` modules  
**Tested in:** `tests/test_integration.py`

---

## F-10 — Nearest-Neighbour Bonds

**Description:** Compute NN distance and coordination number from geometry;
draw bonds in the slab (including bonds whose projections overlap), indicating
in-plane vs. out-of-plane bonds.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F10-AC1 | FCC NN distance = a/√2 (computed, not hardcoded) to 1e-10. |
| F10-AC2 | FCC CN = 12 (computed from geometry). |
| F10-AC3 | BCC NN distance = a√3/2 (computed) to 1e-10. |
| F10-AC4 | BCC CN = 8 (computed from geometry). |
| F10-AC5 | Bonds are drawn between all NN pairs within the slab. |
| F10-AC6 | In-plane bonds (both atoms in the same layer) are visually distinct from out-of-plane. |
| F10-AC7 | No bond is drawn between atoms farther than d_NN + ε. |

**Implemented in:** `core/neighbors.py` + `viz/plot2d.py`  
**Tested in:** `tests/test_neighbors.py`

---

## F-11 — Export

**Description:** Export plots as PNG/SVG, atom data as CSV, 3D coordinates as
XYZ, and computed metrics as JSON.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F11-AC1 | PNG file is a valid image (can be re-opened with matplotlib). |
| F11-AC2 | SVG file is valid XML. |
| F11-AC3 | CSV columns: index, x, y, z, layer, proj_x, proj_y. |
| F11-AC4 | CSV values match the in-memory arrays to machine precision. |
| F11-AC5 | XYZ file follows standard format; re-imported positions match to 1e-10. |
| F11-AC6 | JSON contains keys: d_hkl, planar_density, stacking, nn_distance, cn. |
| F11-AC7 | JSON values match core/ outputs to 1e-10. |

**Implemented in:** `export/exporters.py`  
**Tested in:** `tests/test_integration.py` (Phase 3 gate check)

---

## F-12 — Bicrystal Builder

**Description:** Construct a bicrystal from two grains separated by an (hkl)
boundary plane; analyse the CSL and report Sigma, misorientation, and boundary type.

### F-12a — Construction

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F12-AC1 | Grain 1 atoms satisfy d ≤ 0; Grain 2 atoms satisfy d > 0. |
| F12-AC2 | R^T R = I to 1e-12; det R = +1 to 1e-12. |
| F12-AC3 | CSL preset Sigma5-[100] at 36.87° gives Sigma = 5. |
| F12-AC4 | CSL preset Sigma3-[111] at 60° gives Sigma = 3. |
| F12-AC5 | Rotation by 0° gives a perfect single crystal; removed_count = 0. |
| F12-AC6 | Rotation by 90° about [100] (cubic symmetry) gives a perfect crystal. |
| F12-AC7 | Removal threshold is derived from d_NN, not hardcoded. |
| F12-AC8 | Boundary is classified as tilt/twist/mixed correctly for test cases. |
| F12-AC9 | Sigma is invariant under equivalent axis/angle representations. |

### F-12b — Views

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F12-AC10 | Plan view projects N layers of each grain onto boundary plane using the shared projection engine. |
| F12-AC11 | Each grain's 2D cell and the coincidence cell are drawn correctly. |
| F12-AC12 | Coincidence-cell area = Sigma × A_2D (reference grain). |
| F12-AC13 | Side view (tilt) shows both grains and the boundary line. |

### F-12c — Analysis

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F12-AC14 | Sigma reported correctly for all test cases. |
| F12-AC15 | Misorientation angle reported to within 0.01°. |
| F12-AC16 | Boundary plane indices in both grains' frames are reported. |
| F12-AC17 | Fraction of coincident sites = 1/Sigma (reported to 4 decimal places). |

**Implemented in:** `core/bicrystal.py` + `viz/plot_bicrystal.py`  
**Tested in:** `tests/test_bicrystal.py`

---

## F-13 — Interactive GUI

**Description:** A PySide6 application providing a unified interface for all
other features.

**Acceptance Criteria:**

| ID | Criterion |
|----|-----------|
| F13-AC1 | Structure selector (FCC/BCC) updates plots without restart. |
| F13-AC2 | (hkl) inputs accept integers; invalid inputs show an error message. |
| F13-AC3 | Slab thickness control updates the plot in real time. |
| F13-AC4 | Layer-colour, bond, and unit-cell-overlay toggles function correctly. |
| F13-AC5 | Multi-plane comparison panel shows ≥ 3 planes simultaneously. |
| F13-AC6 | Bicrystal panel: axis, angle, CSL preset, boundary plane, N layers, plan/side view. |
| F13-AC7 | All export buttons produce correct files (verified by re-import). |
| F13-AC8 | Displayed metrics (d, density, stacking, Sigma) match core/ outputs exactly. |
| F13-AC9 | GUI contains no mathematical computations (all calls delegate to core/). |
| F13-AC10 | Application closes cleanly with no unhandled exceptions. |

**Implemented in:** `gui/app.py`  
**Tested in:** `tests/test_integration.py` (Phase 4 scripted cases)
