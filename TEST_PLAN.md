# TEST_PLAN — Validation Strategy and Reference Values

> This document is the authoritative test specification. It is implemented
> as automated pytest tests in `tests/`. Mathematical derivations are in
> MATH_SPEC.md; reference values are in SOURCE_OF_TRUTH.md §6.
>
> **Process rule:** Run the full test suite at the end of every phase and
> report results. If a test fails, diagnose the mathematical cause and fix
> the implementation. Never relax tolerances or edit expected values to
> match output without first proving the expected value wrong from MATH_SPEC.md.

---

## §1 Tolerance Policy

| Category | Tolerance | Rationale |
|----------|-----------|-----------|
| Float comparisons in tests | `1e-12` | Near machine epsilon for double |
| Exact (layer index, Sigma, CN) | Exact integer equality | These must be exact |
| Planar density reference values | `1e-10` | Allows accumulation in area computation |
| Misorientation angle | `0.01°` | Inverse-trig precision |

**Tolerance creep is forbidden.** If a test fails at `1e-12`, the code is wrong,
not the tolerance.

---

## §2 Unit Tests — Phase 1

### §2.1 `test_plane.py`

```
TP-01  validate_and_reduce(0, 0, 0)         → ValueError
TP-02  validate_and_reduce(2, 4, 2)         → (1, 2, 1)
TP-03  validate_and_reduce(-1, -1, -1)      → (-1, -1, -1)
TP-04  unit_normal(1, 1, 1)                 → (1,1,1)/√3  [|n|=1 to 1e-12]
TP-05  unit_normal(h,k,l) · (h,k,l) > 0   [parallel, same sign]
TP-06  unit_normal(-h,-k,-l) = -unit_normal(h,k,l)
TP-07  |unit_normal(h,k,l)| = 1.0 for 100 random (h,k,l)  [hypothesis]
```

### §2.2 `test_projection.py`

```
TP-08  orthonormal_basis(1,1,1): |e1|=|e2|=1, e1·e2=0, e1·n=0, e2·n=0
TP-09  e1 × e2 = n̂  (right-handedness)  for (1,1,1),(1,0,0),(1,1,0)
TP-10  For 100 random (h,k,l): all orthonormality conditions hold [hypothesis]
TP-11  project_2d: for each atom, its 3D position reconstructed as
       proj_x*e1 + proj_y*e2 + d*n̂ equals the original to 1e-12
TP-12  signed_distance_exact: layer key = 2*(h*x+k*y+l*z) is exact integer
       for all FCC and BCC basis atoms
```

### §2.3 `test_layers.py`

```
TP-13  assign_layers(FCC, 1,1,1, 3×3×3 supercell): layer keys are
       contiguous integers; no float comparison used
TP-14  FCC (111) 2×2×2 cell: 6 distinct layers (d=0, 1/√3·a/3 steps)
TP-15  Filter by n_layers=3: exactly 3 distinct layers returned per side
TP-16  Filter by distance range [0, 0.5]: all returned atoms satisfy 0 ≤ d ≤ 0.5a
TP-17  stacking_sequence(FCC, 1,1,1) = ['A','B','C'] (period-3 repeat)
TP-18  stacking_sequence(FCC, 1,0,0) = ['A','B'] (period-2 repeat)
TP-19  stacking_sequence(FCC, 1,1,0) = ['A','B']
TP-20  stacking_sequence(BCC, 1,1,0) = ['A','B']
TP-21  stacking_sequence(BCC, 1,0,0) = ['A','B']
TP-22  stacking_sequence(BCC, 1,1,1) = ['A','B','C']
TP-23  (-h,-k,-l) gives same set of distinct layer distances as (h,k,l)
```

### §2.4 `test_cell2d.py`

```
TP-24  find_inplane_vectors(structure, h,k,l): all returned vectors satisfy
       hw_1 + kw_2 + lw_3 = 0 (integer perpendicularity, w_i = 2u_i);
       parity filter per MATH_SPEC §7.1:
         SC:  all w_i even
         FCC: (w_1+w_2+w_3) % 2 == 0  [sum-even, NOT same-parity]
         BCC: w_1%2 == w_2%2 == w_3%2  [all same parity]
TP-25  Cell vectors are linearly independent: (v_1 x v_2) != (0,0,0)
       tested in exact Fraction arithmetic;
       collinear multiples (2*v_1, (3/2)*v_1, etc.) are correctly rejected
TP-26  gauss_reduce returns |v1|² ≤ |v2|² and v1·v2 / |v1|² ∈ (-0.5, 0.5]
       (Minkowski-reduced pair, tested with exact Fraction arithmetic)
TP-27  cell_area(FCC, 1,1,1) = √3/4 ≈ 0.4330 a²
       [derived: rho=4/√3, so A_2D = 1/rho = √3/4; tested to tol 1e-10]
TP-28  For 50 random (h,k,l) over FCC and BCC, verify rho_hkl × A_2D = 1 via
       BOTH exact squared Fraction check AND float check:
         A_sq  = Fraction( (v1·v1)*(v2·v2) - (v1·v2)^2 )  # exact rational
         d_sq  = Fraction( Delta_d_int**2, 4*(h**2+k**2+l**2) )  # exact rational
         rho_sq = Fraction(N_conv**2) * d_sq
         assert rho_sq * A_sq == Fraction(1, 1)              # exact equality
         assert math.isclose(rho * A_2D, 1.0, rel_tol=1e-12) # float sanity
TP-29  Equivalent planes (permutations under cubic symmetry) give same area
```

### §2.5 `test_density.py`

```
TP-30  formula_spacing(1,0,0) = 1.0 (a=1)
TP-31  formula_spacing(1,1,0) = 1/√2
TP-32  formula_spacing(1,1,1) = 1/√3
TP-33  true_spacing('fcc', 1,1,1) = 1/√3            [no extinction: all odd]
TP-34  true_spacing('fcc', 1,1,0) = 1/(2√2)         [extinction: mixed parity]
TP-35  true_spacing('fcc', 2,0,0) = 1/2              [no extinction: all even]
TP-36  true_spacing('fcc', 1,0,0) = 1/(2·1) = 0.5?
        → Actually (1,0,0) has mixed parity → extinction → true = a/(2√1) = 0.5
        → formula = a/√1 = 1.0 → true = formula/2 = 0.5  ✓
TP-37  true_spacing('bcc', 1,1,0) = 1/√2            [h+k+l=2 even, no extinction]
TP-38  true_spacing('bcc', 1,0,0) = 1/2             [h+k+l=1 odd, extinction]
TP-39  true_spacing('bcc', 1,1,1) = 1/(2√3)         [h+k+l=3 odd, extinction]
TP-40  planar_density('fcc', 1,1,1) = 4/√3 ≈ 2.3094  [tol 1e-10]
TP-41  planar_density('fcc', 1,0,0) = 2.0
TP-42  planar_density('fcc', 1,1,0) = √2 ≈ 1.4142
TP-43  planar_density('bcc', 1,1,0) = √2 ≈ 1.4142
TP-44  planar_density('bcc', 1,0,0) = 1.0
TP-45  planar_density('bcc', 1,1,1) = 1/√3 ≈ 0.5774
TP-46  planar_density(-h,-k,-l) = planar_density(h,k,l)  for all test cases
TP-47  Equivalent planes (hkl) family give same density [at least: (100),(010),(001)]
```

### §2.6 `test_neighbors.py`

```
TP-48  nearest_neighbour('fcc')[0] = 1/√2 ≈ 0.7071  [tol 1e-10]
TP-49  nearest_neighbour('fcc')[1] = 12              [exact integer]
TP-50  nearest_neighbour('bcc')[0] = √3/2 ≈ 0.8660  [tol 1e-10]
TP-51  nearest_neighbour('bcc')[1] = 8               [exact integer]
TP-52  d_NN computed from geometry, not from a hardcoded formula
        (verified by running the same code on a modified basis and checking it adapts)
```

---

## §3 Unit Tests — Phase 2 (Bicrystal)

```
TP-53  rotation_matrix([0,0,1], 0°)  → identity matrix to 1e-12
TP-54  rotation_matrix([1,0,0], 90°): R^T R = I, det R = +1  [tol 1e-12]
TP-55  csl_rotation('sigma5_001'): R^T R = I, det R = +1  [exact Fraction]
       Matrix = (1/5)[[4,-3,0],[3,4,0],[0,0,5]]; axis is [001]
       csl_rotation('sigma5_100'): matrix = (1/5)[[5,0,0],[0,4,-3],[0,3,4]]; axis [100]
TP-56  csl_rotation('sigma3_111'): R^T R = I, det R = +1  [exact Fraction]
       Matrix = (1/3)[[2,-1,2],[2,2,-1],[-1,2,2]]; Tr=2 → θ=60°; axis (1,1,1) ✓
TP-57  compute_sigma(csl_rotation('sigma5_001')) = 5     [exact integer]
       compute_sigma(csl_rotation('sigma5_100')) = 5     [exact integer]
TP-58  compute_sigma(csl_rotation('sigma3_111')) = 3     [exact integer]
TP-59  build_bicrystal(0° rotation): all atoms in one grain, removed_count=0
TP-60  build_bicrystal(90° about [100]): result identical to single crystal
TP-61  Grain 1 atoms all have d ≤ 0; Grain 2 atoms all have d > 0
       [requires centered supercell per MATH_SPEC §1.3]
TP-62  Removal threshold = 0.5 × d_NN; no removed atom has d > 0.5·d_NN to any NN
TP-63  Boundary classification: axis⊥n → 'tilt'; axis∥n → 'twist'
TP-64  Sigma invariant: same axis/angle via different representations → same Sigma
TP-65  Coincidence-cell area in plan view = Sigma × A_2D(reference grain)
TP-66  Fraction of coincident sites in plan view = 1/Sigma
TP-67  Misorientation angle extracted from R matches input angle to 0.01°
```

---

## §4 Integration Tests — Phase 3 (Export Round-Trip)

```
TP-68  Export CSV → reimport → positions match to 1e-10, layer indices exact
TP-69  Export XYZ → reimport → 3D positions match to 1e-10
TP-70  Export JSON → reimport → all numeric values match to 1e-10
TP-71  Export PNG → file exists and size > 0 bytes; re-openable by PIL/matplotlib
TP-72  Export SVG → valid XML (xml.etree.ElementTree parses without error)
```

---

## §5 Integration Tests — Phase 4 (GUI)

```
TP-73  GUI: set structure=FCC, hkl=(1,1,1) → displayed density matches
       planar_density('fcc',1,1,1) exactly
TP-74  GUI: set structure=BCC, hkl=(1,1,0) → displayed stacking matches
       stacking_sequence('bcc',1,1,0)
TP-75  GUI: set CSL preset Sigma5 → displayed Sigma = 5
TP-76  GUI: export CSV from GUI matches direct core/ CSV export bit-for-bit
TP-77  GUI: changing hkl refreshes the plot without crashing
```

---

## §6 Property-Based Tests (Hypothesis)

These run for ≥ 100 random inputs unless stated otherwise.

```
TP-PB-1  Random (h,k,l) ≠ (0,0,0): |unit_normal| = 1  [tol 1e-12]
TP-PB-2  Random (h,k,l): e1·e2=0, |e1|=|e2|=1, e1×e2=n̂  [tol 1e-12]
TP-PB-3  Random (h,k,l): all projected atoms satisfy n̂·r_perp_3D = 0
         where r_perp_3D = (r·e1)*e1 + (r·e2)*e2  (3D in-plane component)
         [MATH_SPEC §14 Invariant 2; tol 1e-12]
TP-PB-4  Random (h,k,l): layer assignment is exact (integer layer key for all atoms)
TP-PB-5  Random (h,k,l) over FCC and BCC, verify rho × A_2D = 1 via exact
         squared Fraction arithmetic (MATH_SPEC §7.4):
           A_sq   = Fraction( (v1·v1)*(v2·v2) - (v1·v2)**2 )
           d_sq   = Fraction(Delta_d_int**2, 4*(h**2+k**2+l**2))
           rho_sq = Fraction(N_conv**2) * d_sq
           assert rho_sq * A_sq == Fraction(1, 1)              # exact
           assert math.isclose(rho * A_2D, 1.0, rel_tol=1e-12) # float sanity
TP-PB-6  Random axis + angle: R^T R = I and det R = +1  [tol 1e-12]
TP-PB-7  Random (h,k,l): (hkl) and (-h,-k,-l) give same d_true and density.
         Stacking assertion is basis-independent (MATH_SPEC §14 Invariant 5):
           P_pos = stacking_period_length(h,k,l)
           P_neg = stacking_period_length(-h,-k,-l)
           assert P_pos == P_neg                               # exact integer
           seq_pos = canonical_stacking_string(h,k,l)         # e.g. 'ABC'
           seq_neg = canonical_stacking_string(-h,-k,-l)      # ascending d
           assert seq_neg == canonical_relabel(reverse(seq_pos))
         [Do NOT assert 'identical set of reduced offsets'; negating n_hat
          applies a unimodular transformation that changes raw (s1,s2) pairs]
```

---

## §7 Phase Gate Summary

| Phase | Tests that must pass | Gate command |
|-------|---------------------|-------------|
| 1 | TP-01 → TP-52 + TP-PB-1 → TP-PB-7 | `pytest tests/ -k "not bicrystal and not integration"` |
| 2 | TP-53 → TP-67 (bicrystal unit tests) | `pytest tests/test_bicrystal.py` |
| 3 | TP-68 → TP-72 (export round-trip) | `pytest tests/test_integration.py -k "export"` |
| 4 | TP-73 → TP-77 (GUI scripted cases) | `pytest tests/test_integration.py -k "gui"` |
| 5 | Full suite | `pytest tests/ -v` |

---

## §8 Traceability Matrix (Phase 5)

| Test ID(s) | Feature | MATH_SPEC section | FEATURES.md acceptance criteria |
|------------|---------|-------------------|--------------------------------|
| TP-01–07 | F-01 | §3.1–3.3 | F01-AC1–6 |
| TP-08–12 | F-02 | §6 | F02-AC1–7 |
| TP-13–23 | F-03, F-05, F-08 | §4, §9 | F03-AC1–5, F05-AC1–4, F08-AC1–8 |
| TP-24–29 | F-07 | §7 | F07-AC1–6 |
| TP-30–47 | F-06 | §5, §8 | F06-AC1–7 |
| TP-48–52 | F-10 | §10 | F10-AC1–4 |
| TP-53–67 | F-12 | §11–13 | F12-AC1–17 |
| TP-68–72 | F-11 | — | F11-AC1–7 |
| TP-73–77 | F-13 | — | F13-AC1–10 |
| TP-PB-1–7 | F-01,02,03,06 | §3,6,8,14 | Various |

---

## §9 Reference Values for Verification

All in units of `a = 1`.

### Planar Density

| Structure | (hkl) | Expected | Exact form |
|-----------|-------|----------|------------|
| FCC | (111) | 2.30940107676… | 4/√3 |
| FCC | (100) | 2.0 | 2 |
| FCC | (110) | 1.41421356237… | √2 |
| BCC | (110) | 1.41421356237… | √2 |
| BCC | (100) | 1.0 | 1 |
| BCC | (111) | 0.57735026919… | 1/√3 |

### True Interplanar Spacing

| Structure | (hkl) | d_formula | d_true | Extinct? |
|-----------|-------|-----------|--------|----------|
| FCC | (111) | 0.57735… | 0.57735… | No (all odd) |
| FCC | (100) | 1.0 | 0.5 | Yes (mixed) |
| FCC | (110) | 0.70711… | 0.35355… | Yes (mixed) |
| FCC | (200) | 0.5 | 0.5 | No (all even) |
| BCC | (110) | 0.70711… | 0.70711… | No (sum=2 even) |
| BCC | (100) | 1.0 | 0.5 | Yes (sum=1 odd) |
| BCC | (111) | 0.57735… | 0.28868… | Yes (sum=3 odd) |

### Nearest-Neighbour

| Structure | d_NN | CN |
|-----------|------|-----|
| FCC | 0.70711… (= 1/√2) | 12 |
| BCC | 0.86603… (= √3/2) | 8 |

### Bicrystal Sigma

| Preset | Axis | Angle | Σ | Exact matrix |
|--------|------|-------|---|---------------|
| sigma5_001 | [001] | 36.87° | 5 | (1/5)[[4,-3,0],[3,4,0],[0,0,5]] |
| sigma5_100 | [100] | 36.87° | 5 | (1/5)[[5,0,0],[0,4,-3],[0,3,4]] |
| sigma3_111 | [111] | 60° | 3 | (1/3)[[2,-1,2],[2,2,-1],[-1,2,2]] |

> **Note on sigma3_111 matrix:** The previously documented matrix
> `(1/3)[[1,-2,2],[2,1,-2],[-2,2,1]]` has Tr=1 implying θ=90° about [011] —
> this is incorrect. The matrix above has Tr=2 (θ=60°) and fixed axis [111]. ✓
