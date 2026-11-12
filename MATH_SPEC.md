# MATH_SPEC — Mathematical Specification

> **First-principles derivations precede every formula.**
> This document is the mathematical ground truth for all implementation decisions.
> No formula appears in code without a reference to a section here.

---

## §1 Bravais Lattice and Conventional Cubic Cell

### §1.1 Primitive vs. Conventional Cell

A Bravais lattice is the infinite set of translation vectors

$$\mathbf{R} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + n_3 \mathbf{a}_3, \quad n_i \in \mathbb{Z}$$

For a **cubic** system with parameter $a$, we work with the **conventional** cubic
cell, which has **orthogonal** lattice vectors:

$$\mathbf{a}_1 = a\hat{x}, \quad \mathbf{a}_2 = a\hat{y}, \quad \mathbf{a}_3 = a\hat{z}$$

The conventional cell is not primitive for FCC and BCC; it contains multiple
atoms per cell (the **basis**).

### §1.2 Bases

**FCC** (face-centred cubic) — 4 atoms per conventional cell:

$$\mathbf{\tau}_0 = (0, 0, 0)a, \quad
\mathbf{\tau}_1 = \tfrac{a}{2}(\hat{x}+\hat{y}), \quad
\mathbf{\tau}_2 = \tfrac{a}{2}(\hat{x}+\hat{z}), \quad
\mathbf{\tau}_3 = \tfrac{a}{2}(\hat{y}+\hat{z})$$

In fractional coordinates $\mathbf{r} = (u, v, w)$ with $u,v,w \in [0,1)$:

$$\{(0,0,0),\ (1/2,1/2,0),\ (1/2,0,1/2),\ (0,1/2,1/2)\}$$

**BCC** (body-centred cubic) — 2 atoms per conventional cell:

$$\mathbf{\tau}_0 = (0, 0, 0)a, \quad \mathbf{\tau}_1 = \tfrac{a}{2}(\hat{x}+\hat{y}+\hat{z})$$

In fractional coordinates:

$$\{(0,0,0),\ (1/2,1/2,1/2)\}$$

### §1.3 Supercell Generation

To generate a symmetric slab centred on the origin, enumerate all conventional
cells in the range:

$$\mathbf{r}_{(n_1,n_2,n_3,b)} = (n_1 + \tau_b^x)\hat{x} + (n_2 + \tau_b^y)\hat{y} + (n_3 + \tau_b^z)\hat{z}$$

for $n_i \in \{-N_i, -N_i+1, \ldots, N_i\}$ and $b$ ranging over basis atoms. This
centres the supercell on the origin so that both $d \le 0$ and $d > 0$ half-spaces
are populated equally for any $(hkl)$, and all projection views are symmetric.
All positions are stored in units of $a$ (i.e., $a=1$ internally).

### §1.4 Aperture Cutoff

To produce a symmetric circular field of view, apply one of the following
cutoffs **after** generating the full centred supercell:

- **Spherical aperture:** retain atoms with $|\mathbf{r}| \le R_\text{max}$.
- **Cylindrical aperture:** retain atoms with $|\mathbf{r}_\perp^{3D}| \le R_\text{max}$
  **and** $|d| \le d_\text{max}$, where $\mathbf{r}_\perp^{3D} = \mathbf{r} - d\hat{n}$.

When rendering projected labels (e.g. fractional-coordinate annotations), two
atoms that project to the **exact same 2D coordinate** $(\mathbf{r}\cdot\mathbf{e}_1,
\mathbf{r}\cdot\mathbf{e}_2)$ are deduplicated by retaining only the atom with the
smallest $|d|$ (i.e., the one closest to the boundary plane).

---

## §2 Reciprocal Lattice

For a lattice with primitive vectors $\{\mathbf{a}_i\}$, the reciprocal vectors $\{\mathbf{b}_i\}$
satisfy:

$$\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi \delta_{ij}$$

For the **conventional cubic** cell ($\mathbf{a}_i = a\hat{e}_i$):

$$\mathbf{b}_1 = \frac{2\pi}{a}\hat{x}, \quad
\mathbf{b}_2 = \frac{2\pi}{a}\hat{y}, \quad
\mathbf{b}_3 = \frac{2\pi}{a}\hat{z}$$

The reciprocal lattice vectors $\{\mathbf{b}_i\}$ are stored internally, and all
plane-normal calculations use them. This keeps the code valid for general lattices.

---

## §3 Plane (hkl) and Its Normal

### §3.1 Definition

A **Miller plane** $(hkl)$ is indexed relative to the **reciprocal lattice**. Its
normal (un-normalised) is the reciprocal-space vector:

$$\mathbf{G}_{hkl} = h\mathbf{b}_1 + k\mathbf{b}_2 + l\mathbf{b}_3$$

For the cubic case this simplifies to:

$$\mathbf{G}_{hkl} = \frac{2\pi}{a}(h, k, l)$$

The direction of $\mathbf{G}_{hkl}$ is perpendicular to the family of $(hkl)$ planes —
this is a consequence of the definition of the reciprocal lattice.

### §3.2 Unit Normal

$$\hat{n} = \frac{\mathbf{G}_{hkl}}{|\mathbf{G}_{hkl}|} = \frac{(h, k, l)}{\sqrt{h^2+k^2+l^2}}$$

(The $2\pi/a$ factors cancel.) This is the unit vector used for all subsequent
projections.

### §3.3 Validation

- $(h, k, l) = (0, 0, 0)$ is **rejected** (no well-defined plane).
- $(h, k, l)$ is **reduced** by $g = \gcd(|h|, |k|, |l|)$ to canonical form.
- Negatives are allowed; $(\bar{h}\bar{k}\bar{l})$ is the antipodal plane.

---

## §4 Signed Distance

The signed distance from the reference plane (the plane through the origin with
normal $\hat{n}$) to an atom at position $\mathbf{r}$ is:

$$d = \hat{n} \cdot \mathbf{r}$$

For the cubic case with $\mathbf{r} = (x, y, z)$ (in units of $a$):

$$d = \frac{hx + ky + lz}{\sqrt{h^2+k^2+l^2}}$$

**Exact computation:** Since we work in units of $a$, coordinates are of the form
$n/2$ for FCC (with $n \in \mathbb{Z}$). To use exact arithmetic, multiply through by
$2$ to get integer coordinates and compute the integer dot product
$2d' = 2(hx + ky + lz) \in \mathbb{Z}$. The layer index is then $2d' \in \mathbb{Z}$.
We store layer indices as integers (the $2\times$ convention) or as `fractions.Fraction`.

---

## §5 True Interplanar Spacing vs. Formula Spacing

### §5.1 Formula Spacing

The standard crystallographic formula gives the spacing between adjacent
*lattice* planes (not necessarily all occupied by atoms):

$$d_{hkl}^{\text{formula}} = \frac{a}{\sqrt{h^2+k^2+l^2}} \quad \text{(cubic)}$$

### §5.2 Structure Factor and Extinction

The structure factor for a reflection $(hkl)$ is:

$$F_{hkl} = \sum_b f_b \exp\!\left(2\pi i (h\tau_b^x + k\tau_b^y + l\tau_b^z)\right)$$

where the sum is over basis atoms. When $F_{hkl} = 0$, the reflection is
**systematically absent** (extinct) and the atom planes at that lattice spacing
are **interleaved** — producing additional intermediate planes that halve the
true interplanar spacing.

**FCC extinction rule:** $F_{hkl} = 0$ when $h, k, l$ are of **mixed parity**
(some even, some odd). In this case the true atomic-plane spacing is:

$$d_{hkl}^{\text{true}} = \frac{a}{2\sqrt{h^2+k^2+l^2}} \quad \text{(mixed parity, FCC)}$$

**BCC extinction rule:** $F_{hkl} = 0$ when $h + k + l$ is **odd**. In this case:

$$d_{hkl}^{\text{true}} = \frac{a}{2\sqrt{h^2+k^2+l^2}} \quad \text{(h+k+l odd, BCC)}$$

### §5.3 Derivation from Atomic Layers

The correct approach (and the one used in code) is to **measure** $d^{\text{true}}$
directly from the distinct $d$-values of the atomic layers:

1. Compute all signed distances $d_i = \hat{n} \cdot \mathbf{r}_i$.
2. Sort and deduplicate (using exact arithmetic).
3. The true interplanar spacing is $\min_i(d_{i+1} - d_i)$ over the periodic
   set of distinct layers.

This measured value must equal $d^{\text{formula}}$ for all-same-parity FCC and
even-sum BCC, and $d^{\text{formula}}/2$ for the extinction cases.

**Test:** Assert that the measured spacing matches the predicted formula (with
extinction correction) to machine precision for all reference planes.

---

## §6 Orthonormal In-Plane Basis

### §6.1 Construction

Given the unit normal $\hat{n}$, we need two orthonormal vectors $\mathbf{e}_1$,
$\mathbf{e}_2$ spanning the plane such that $(\mathbf{e}_1, \mathbf{e}_2, \hat{n})$
is **right-handed**.

**Step 1:** Find a vector not parallel to $\hat{n}$. We use the axis of the
smallest component of $\hat{n}$:

$$\mathbf{t} = \hat{x} \quad \text{if } |n_x| \le |n_y|, |n_z|, \quad
\text{else} \quad \mathbf{t} = \hat{y} \quad \text{or} \quad \hat{z}$$

**Step 2:** Gram-Schmidt:

$$\mathbf{e}_1 = \frac{\mathbf{t} - (\mathbf{t}\cdot\hat{n})\hat{n}}{|\mathbf{t} - (\mathbf{t}\cdot\hat{n})\hat{n}|}$$

**Step 3:** Right-hand rule:

$$\mathbf{e}_2 = \hat{n} \times \mathbf{e}_1$$

**Verification:** Check $|\mathbf{e}_1| = |\mathbf{e}_2| = 1$, $\mathbf{e}_1\cdot\mathbf{e}_2 = 0$,
$\mathbf{e}_1\cdot\hat{n} = \mathbf{e}_2\cdot\hat{n} = 0$,
$\mathbf{e}_1 \times \mathbf{e}_2 = \hat{n}$ (right-handedness).

### §6.2 Orthogonal Projection to 2D

The 2D projected position of an atom $\mathbf{r}$ onto the $(hkl)$ plane is:

$$\mathbf{r}_\perp = (\mathbf{r}\cdot\mathbf{e}_1,\ \mathbf{r}\cdot\mathbf{e}_2)$$

The residual $\mathbf{r} - (\mathbf{r}\cdot\mathbf{e}_1)\mathbf{e}_1 - (\mathbf{r}\cdot\mathbf{e}_2)\mathbf{e}_2$
equals $(\hat{n}\cdot\mathbf{r})\hat{n}$, confirming the atom lies at signed distance $d$
from the reference plane.

---

## §7 2D In-Plane Lattice and Primitive Cell

### §7.1 True 2D Bravais In-Plane Lattice Vectors

An in-plane lattice vector must be a **translation of the 3D Bravais lattice**
that is **perpendicular to $\hat{n}$**. For FCC and BCC, the primitive 3D Bravais
translations include half-integer steps that the simple-cubic integer grid misses.
Using only $u_i \in \mathbb{Z}$ (the simple-cubic grid) finds the wrong (too large)
2D cell for FCC and BCC, breaking §7.4 and §9.2.

The correct candidate set for each structure is expressed using **doubled integer coordinates** $w_i = 2u_i \in \mathbb{Z}$, which avoids fractions in the parity tests:

| Structure | Allowed $\mathbf{u}$ | Doubled-coord rule on $\mathbf{w} = (w_1,w_2,w_3) = (2u_1,2u_2,2u_3)$ |
|-----------|---------------------|-------------------------------------------------------------------------|
| **SC** | $u_i \in \mathbb{Z}$ | $w_1, w_2, w_3$ **all even**: $w_i \equiv 0 \pmod{2}$ |
| **FCC** | $u_i \in \tfrac{1}{2}\mathbb{Z}$ | $w_1 + w_2 + w_3 \equiv 0 \pmod{2}$ (sum of doubled coords is **even** — either all even, or exactly two odd and one even) |
| **BCC** | $u_i \in \tfrac{1}{2}\mathbb{Z}$ | $w_1 \equiv w_2 \equiv w_3 \pmod{2}$ (all three doubled coords have the **same parity** — all even or all odd) |

> **Why these rules?**
> - **SC:** The primitive SC cell has translations $n_1\mathbf{a}_1+n_2\mathbf{a}_2+n_3\mathbf{a}_3$ with $n_i \in \mathbb{Z}$, so $w_i = 2n_i$ are all even.
> - **FCC:** A valid FCC Bravais translation is any sum of a lattice vector plus any face-centering offset. The four FCC basis vectors have $\mathbf{w}$ values $(0,0,0)$, $(1,1,0)$, $(1,0,1)$, $(0,1,1)$ — all with **even sum**. Any integer linear combination of these also has even $w_1+w_2+w_3$. For example, $\mathbf{u}=(1/2,1/2,0)$ gives $\mathbf{w}=(1,1,0)$, sum $= 2$ (even) ✓; whereas $\mathbf{u}=(1/2,1/2,1/2)$ gives $\mathbf{w}=(1,1,1)$, sum $= 3$ (odd) ✗ — this is a BCC body-centre, not an FCC translation.
> - **BCC:** The BCC primitive translation is $(1/2,1/2,1/2)$, giving $\mathbf{w}=(1,1,1)$ — all **same parity** (all odd). The general BCC translation $(n_1 + m/2, n_2 + m/2, n_3 + m/2)$ with $m \in \{0,1\}$ gives $w_i$ all even (m=0) or all odd (m=1).

In practice, enumerate candidates with $|w_i| \le 2B$ for some bounding radius $B$, apply the structure's parity filter in integer arithmetic, then filter by perpendicularity:

$$hu_1 + ku_2 + lu_3 = 0 \iff \frac{h w_1 + k w_2 + l w_3}{2} = 0 \iff hw_1 + kw_2 + lw_3 = 0$$

All arithmetic uses **exact integers** (or `fractions.Fraction` for $u_i = w_i/2$), so the perpendicularity test is exact.

> **Result:** Because each primitive 3D cell of FCC (or BCC) contains exactly
> 1 atom, the 2D primitive layer cell also contains exactly 1 atom, so
> $\rho_{hkl} \times A_{2D} = 1$ for all $(hkl)$. The consistency check in §7.4
> therefore collapses to the integer 1.

### §7.2 Gauss/Lagrange Reduction

**Selecting the initial pair.** From the enumerated candidate vectors (sorted by
squared length $|\mathbf{u}|^2 = \mathbf{u} \cdot \mathbf{u}$, computed exactly as
`Fraction`), pick $\mathbf{v}_1$ as the shortest non-zero candidate, then pick
$\mathbf{v}_2$ as the next shortest candidate that is **linearly independent** of
$\mathbf{v}_1$. Linear independence is tested in exact `Fraction` arithmetic via the
3D cross product:

$$\mathbf{v}_1 \times \mathbf{v}_2 \ne (0, 0, 0)$$

> **Why not just check $\mathbf{v}_2 \ne \pm\mathbf{v}_1$?** Non-primitive multiples
> such as $2\mathbf{v}_1$ or $\tfrac{3}{2}\mathbf{v}_1$ are parallel to $\mathbf{v}_1$
> and will appear in the candidate list whenever the bounding radius is large
> enough. The cross-product test rejects all such collinear vectors exactly.

**Gauss/Lagrange reduction.** Given the two selected independent vectors, apply
reduction to find the provably shortest (Minkowski-reduced) pair. All arithmetic
is in exact `Fraction` coordinates:

```
Gauss reduction (exact):
  while |v2|² < |v1|²: swap v1, v2          # compare squared lengths exactly
  n = floor_round(v1·v2 / |v1|²)             # nearest integer (exact Fraction → int)
  v2 = v2 - n*v1
  repeat until no swap and n == 0
```

The squared lengths $|\mathbf{v}|^2$ and dot products $\mathbf{v}_1 \cdot \mathbf{v}_2$ are
computed exactly as `Fraction` values (since $v_i \in \tfrac{1}{2}\mathbb{Z}$).

The result is the pair of shortest independent in-plane Bravais translations,
defining the **2D primitive cell**.

### §7.3 2D Cell Area

$$A_{2D} = |\mathbf{v}_1 \times \mathbf{v}_2| \quad \text{(magnitude of 3D cross product)}$$

Equivalently, in 2D projected coordinates:

$$A_{2D} = |v_{1,x}^{2D} v_{2,y}^{2D} - v_{1,y}^{2D} v_{2,x}^{2D}|$$

This is computed from float projections at the render boundary; the exact area
can also be derived from the `Fraction` 3D vectors.

### §7.4 Consistency Check

Since each primitive 3D Bravais cell contains exactly 1 atom (SC: 1, FCC: 1,
BCC: 1 — the conventional cell multi-atom count is exactly cancelled by the
larger conventional-cell volume), the 2D primitive cell derived from §7.1
also contains exactly **1 atom**, so:

$$\rho_{hkl} \times A_{2D} = 1 \quad \text{(exactly, for SC / FCC / BCC)}$$

**Implementation note — exact rational verification via squared quantities.**
Both $\rho_{hkl}$ and $A_{2D}$ contain the irrational factor $1/\sqrt{h^2+k^2+l^2}$,
so `fractions.Fraction` cannot represent them directly, and IEEE-754 float
multiplication is subject to rounding. The exact test is therefore performed on
the **squared** quantities, which are rational:

$$A_{2D}^2 = |\mathbf{v}_1 \times \mathbf{v}_2|^2
           = (\mathbf{v}_1 \cdot \mathbf{v}_1)(\mathbf{v}_2 \cdot \mathbf{v}_2)
             - (\mathbf{v}_1 \cdot \mathbf{v}_2)^2 \in \mathbb{Q}$$

$$\left(d_{hkl}^\text{true}\right)^2
  = \frac{(\Delta d_\text{int})^2}{4(h^2+k^2+l^2)} \in \mathbb{Q}$$

$$\rho_{hkl}^2 = N_\text{conv}^2 \cdot \left(d_{hkl}^\text{true}\right)^2 \in \mathbb{Q}$$

where $\Delta d_\text{int}$ is the integer numerator of the minimum layer spacing
(from §5.3). The test suite asserts **both**:

1. **Exact rational equality** (using `Fraction`):
   $$\rho_{hkl}^2 \times A_{2D}^2 = \texttt{Fraction}(1, 1)$$
2. **Float consistency** (as a secondary sanity check):
   `math.isclose(rho_hkl * A_2D, 1.0, rel_tol=1e-12)`

---

## §8 Planar Density

### §8.1 First-Principles Derivation

The planar density $\rho_{hkl}$ is the number of atoms per unit area in one
atomic plane of type $(hkl)$.

**Derivation from volume per atom:**
A conventional cubic cell of side $a$ has volume $a^3$ and contains $N_\text{conv}$
atoms ($N_\text{conv} = 4$ for FCC, $2$ for BCC). The volume per atom is therefore:

$$\Omega = \frac{a^3}{N_\text{conv}}$$

In a crystal, the total volume is tiled by atomic layers of area $A_{hkl}$ spaced
$d_{hkl}^\text{true}$ apart. The volume per atom assigned to one layer is therefore
$A_{hkl} \times d_{hkl}^\text{true}$ divided by the number of atoms in area $A_{hkl}$.
Since each atom occupies volume $\Omega$ and the layer spacing allocates a column
of height $d_{hkl}^\text{true}$ per atom, we have:

$$\rho_{hkl} \times d_{hkl}^\text{true} = \frac{1}{\Omega} = \frac{N_\text{conv}}{a^3}$$

Solving:

$$\boxed{\rho_{hkl} = \frac{N_\text{conv}\, d_{hkl}^\text{true}}{a^3}}$$

With $a = 1$ internally:

$$\rho_{hkl} = N_\text{conv}\, d_{hkl}^\text{true}$$

This formula is exact and first-principles; it requires only $N_\text{conv}$ and the
true interplanar spacing $d_{hkl}^\text{true}$ from §5.

**Consistency with §7.4:** Since the primitive 3D Bravais cell has volume
$\Omega = a^3/N_\text{conv}$ and contains exactly 1 atom, the 2D primitive layer
cell area is $A_{2D} = \Omega / d_{hkl}^\text{true} = 1 / \rho_{hkl}$,
which gives $\rho_{hkl} \times A_{2D} = 1$ (exactly).

### §8.2 Reference Values (units: $a^{-2}$, with $a = 1$)

| Structure | $(hkl)$ | $d^\text{true}$ | $\rho_{hkl}$ | Exact form |
|-----------|---------|-----------------|--------------|------------|
| FCC | $(111)$ | $1/(2\sqrt{3}) \times 2 = 1/\sqrt{3}$\* | $4/\sqrt{3} \approx 2.309$ | $4/\sqrt{3}$ |
| FCC | $(100)$ | $1/2$ | $2$ | $2$ |
| FCC | $(110)$ | $1/(2\sqrt{2})$ | $\sqrt{2} \approx 1.414$ | $\sqrt{2}$ |
| BCC | $(110)$ | $1/\sqrt{2}$ | $\sqrt{2} \approx 1.414$ | $\sqrt{2}$ |
| BCC | $(100)$ | $1/2$ | $1$ | $1$ |
| BCC | $(111)$ | $1/(2\sqrt{3})$ | $1/\sqrt{3} \approx 0.577$ | $1/\sqrt{3}$ |

\*FCC (111) has all-odd indices so there is **no extinction**; the true spacing
is $a/\sqrt{3}$ (not halved), giving $\rho = 4 \times (1/\sqrt{3}) = 4/\sqrt{3}$.

Displayed as explicit equations:

$$\rho_\text{FCC,111} = \frac{4}{\sqrt{3}} \approx 2.309 \,a^{-2}$$
$$\rho_\text{FCC,100} = 2 \,a^{-2}$$
$$\rho_\text{FCC,110} = \sqrt{2} \approx 1.414 \,a^{-2}$$
$$\rho_\text{BCC,110} = \sqrt{2} \approx 1.414 \,a^{-2}$$
$$\rho_\text{BCC,100} = 1 \,a^{-2}$$
$$\rho_\text{BCC,111} = \frac{1}{\sqrt{3}} \approx 0.577 \,a^{-2}$$

---

## §9 Stacking Sequence

### §9.1 Layer Offsets

For each layer $i$ (sorted by signed distance $d_i$), compute the full 3D
difference vector between a representative atom of layer $i$ and one of layer 0,
working in exact `Fraction` arithmetic throughout:

$$\boldsymbol{\delta}_i^{3D} = \mathbf{r}^{(i)}_\text{rep} - \mathbf{r}^{(0)}_\text{rep}$$

where $\mathbf{r}^{(j)}_\text{rep}$ is any atom in layer $j$ (all atoms in a layer
share the same in-plane offset modulo the 2D lattice). All coordinates are
`fractions.Fraction`. Note that $\boldsymbol{\delta}_i^{3D}$ **includes** the
out-of-plane component $(d_i - d_0)\hat{n}$; this is handled automatically in §9.2.

### §9.2 Exact Fractional Modulo Reduction

Express $\boldsymbol{\delta}_i^{3D}$ as a linear combination of the two 3D in-plane
primitive vectors $\mathbf{v}_1$, $\mathbf{v}_2$ (from §7.1) by solving the **normal
equations** in exact rational arithmetic:

$$\begin{pmatrix} \mathbf{v}_1\cdot\mathbf{v}_1 & \mathbf{v}_1\cdot\mathbf{v}_2 \\
                  \mathbf{v}_2\cdot\mathbf{v}_1 & \mathbf{v}_2\cdot\mathbf{v}_2 \end{pmatrix}
\begin{pmatrix} s_1 \\ s_2 \end{pmatrix}
= \begin{pmatrix} \mathbf{v}_1 \cdot \boldsymbol{\delta}_i^{3D} \\
                  \mathbf{v}_2 \cdot \boldsymbol{\delta}_i^{3D} \end{pmatrix}$$

> **Normal component cancellation.** Although $\boldsymbol{\delta}_i^{3D}$ contains
> the out-of-plane part $(d_i - d_0)\hat{n}$, the right-hand sides
> $\mathbf{v}_k \cdot \boldsymbol{\delta}_i^{3D}$ automatically project it out:
> $$\mathbf{v}_k \cdot \boldsymbol{\delta}_i^{3D}
>   = \mathbf{v}_k \cdot \boldsymbol{\delta}_{i,\perp}^{3D}
>     + (d_i - d_0)\underbrace{(\mathbf{v}_k \cdot \hat{n})}_{=\,0}$$
> because $\mathbf{v}_k$ is in-plane ($\mathbf{v}_k \cdot \hat{n} = 0$ by §7.1).
> The normal equations therefore yield the **exact in-plane** fractional
> coordinates $(s_1, s_2)$ without any need to subtract $(d_i - d_0)\hat{n}$
> beforehand.
>
> *Residual check only:* if one later reconstructs the 3D in-plane offset
> $s_1\mathbf{v}_1 + s_2\mathbf{v}_2$ and wishes to verify it against
> $\boldsymbol{\delta}_i^{3D}$, compare against the in-plane projection:
> $$\boldsymbol{\delta}_{i,\perp}^{3D}
>   = \boldsymbol{\delta}_i^{3D}
>     - \frac{\boldsymbol{\delta}_i^{3D} \cdot \mathbf{G}_{hkl}}
>            {|\mathbf{G}_{hkl}|^2}\,\mathbf{G}_{hkl}$$
> where $\mathbf{G}_{hkl} = (h,k,l)$ (the un-normalised normal, §3.1).

All entries in the Gram matrix are exact `Fraction` values (since
$\mathbf{v}_j \in \tfrac{1}{2}\mathbb{Z}^3$ and $\boldsymbol{\delta}_i \in \tfrac{1}{2}\mathbb{Z}^3$).
Solve by Cramer's rule or Fraction-based Gaussian elimination to obtain
$s_1, s_2 \in \mathbb{Q}$.

Reduce each fractional coordinate into $[0, 1)$ using **floor division**:

$$s_k^\text{mod} = s_k - \lfloor s_k \rfloor \in [0, 1) \cap \mathbb{Q}$$

This is computed exactly with `Fraction` floor: `s_mod = s - Fraction(math.floor(s))`.

> **Why not `round(M^{-1}\delta)`?** For offsets on half-integer cell boundaries
> (e.g., $s = 1/2$ or $s = -1/2$), `round` has implementation-defined tie-breaking
> that can map equivalent offsets to different integers, misidentifying stacking
> labels. The floor-based reduction is **deterministic and exact** for all rational
> inputs.

### §9.3 Label Assignment

Two layers $i$ and $j$ share the same stacking label **if and only if** their
reduced fractional coordinates are **identically equal** as `Fraction` pairs:

$$(s_1^{\text{mod},i},\ s_2^{\text{mod},i}) = (s_1^{\text{mod},j},\ s_2^{\text{mod},j})$$

Labels A, B, C, … are assigned in order of first appearance across layers
$0, 1, 2, \ldots$ (sorted by ascending $d$). The stacking sequence is the
resulting label string; it must repeat with a period equal to the number of
distinct reduced offsets.

---

## §10 Nearest-Neighbour Geometry

### §10.1 Derivation

**FCC:** The 12 nearest neighbours of a corner atom are the 12 face-centre
atoms of the two adjacent cells sharing that corner. Each is at distance:

$$d_\text{NN}^\text{FCC} = \frac{a}{\sqrt{2}}$$

**Derivation:** Face-centre atom $\mathbf{\tau}_1 = (a/2)(1,1,0)$ has distance from origin
$|\mathbf{\tau}_1| = a\sqrt{(1/2)^2+(1/2)^2} = a/\sqrt{2}$. The 12 NN are all
$(\pm 1/2, \pm 1/2, 0)$-type permutations (choosing 2 of 3 axes, 4 sign
combinations each, $\binom{3}{2}\times 4 = 12$). ✓

**BCC:** The 8 nearest neighbours of a corner atom are the body centres of
the 8 cubes sharing that corner:

$$d_\text{NN}^\text{BCC} = \frac{a\sqrt{3}}{2}$$

**Derivation:** Body-centre atom $\mathbf{\tau}_1 = (a/2)(1,1,1)$ has distance
$|\mathbf{\tau}_1| = a\sqrt{3}/2$. All $(\pm 1/2, \pm 1/2, \pm 1/2)$ vectors have the same
length; there are $2^3 = 8$ such vectors. ✓

### §10.2 Computation in Code

We do NOT hardcode these distances. Instead:
1. Build a small supercell (e.g., $3\times 3\times 3$ conventional cells).
2. For each atom, compute distances to all other atoms.
3. Identify the minimum distance $d_\text{min}$ across the supercell (excluding self).
4. Count atoms at distance $d_\text{min}$ (within a tight tolerance; since all
   coordinates are half-integer multiples, distances are $a\sqrt{n}/2$ for integer $n$,
   so exact grouping is possible by comparing $4d^2/a^2 \in \mathbb{Z}$).

---

## §11 Rodrigues' Rotation Formula

For an **active** rotation by angle $\theta$ about unit axis $\hat{k}$, acting on
column vectors $\mathbf{v}' = R\mathbf{v}$:

$$\mathbf{v}' = \mathbf{v}\cos\theta + (\hat{k}\times\mathbf{v})\sin\theta + \hat{k}(\hat{k}\cdot\mathbf{v})(1-\cos\theta)$$

Written as a matrix, the $(i,j)$ entry of $R$ is obtained by evaluating the above
for $\mathbf{v} = \hat{e}_j$ and reading off the $i$-th component:

$$R_{ij} = \cos\theta\,\delta_{ij} + (1-\cos\theta)k_i k_j - \sin\theta\,\varepsilon_{ijk}k_k$$

where $\varepsilon_{ijk}$ is the Levi-Civita symbol (fully antisymmetric, $\varepsilon_{123}=+1$).

> **Sign note:** The term is $-\sin\theta\,\varepsilon_{ijk}k_k$ (minus sign).
> This follows from $(\hat{k} \times \mathbf{v})_i = \varepsilon_{ijk}k_j v_k$,
> so the $j$-th column of $R$ has $i$-th entry
> $\delta_{ij}\cos\theta + (1-\cos\theta)k_i k_j + \varepsilon_{ijk}k_j\sin\theta$.
> Since $\varepsilon_{ijk}k_j = -\varepsilon_{ikj}k_j = -(\hat{k}\times)_{ij}$ when
> written in terms of the Levi-Civita with the third index as the summation index,
> the convention above is internally consistent with the vector formula. In NumPy:
> ```python
> K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])  # skew-sym
> R = np.eye(3)*np.cos(t) + (1-np.cos(t))*np.outer(k,k) + np.sin(t)*K
> ```
> Here `K` encodes $\hat{k}\times$ such that $K\mathbf{v} = \hat{k}\times\mathbf{v}$,
> matching the vector formula with a **plus** sign on `sin(t)*K`.

**Verification:** $R^T R = I$ and $\det R = +1$ must hold (asserted in code, not just tested).

### §11.1 CSL Exact Rotations

For CSL presets, rotation matrices are given as exact rational matrices
(entries in $\mathbb{Q}$), **not** computed from floating-point trigonometry.

**Sigma 5 about [001] (36.87°):**

The rotation by $\theta = \arctan(3/4) \approx 36.87°$ about the $[001]$ axis:

$$R_{\Sigma 5} = \frac{1}{5}\begin{pmatrix} 4 & -3 & 0 \\ 3 & 4 & 0 \\ 0 & 0 & 5 \end{pmatrix}$$

> **Axis clarification:** This matrix has its fixed axis along $[001]$ (the $z$-axis),
> not $[100]$. The equivalent matrix for a rotation about $[100]$ is:
> $\frac{1}{5}\begin{pmatrix}5&0&0\\0&4&-3\\0&3&4\end{pmatrix}$.
> Both generate $\Sigma = 5$. The preset `sigma5_001` uses the first form;
> `sigma5_100` uses the second.

Verification (exact, using `Fraction`):
- $\text{Tr}(R) = (4+4+5)/5 = 13/5$. For a rotation by $\theta$, $\text{Tr}(R) = 1 + 2\cos\theta$. Here $\cos(36.87°) = 4/5$, so $\text{Tr}(R) = 1 + 8/5 = 13/5$. ✓
- $R^T R = I$ (exact). ✓
- $\det R = (4\cdot4-(-3)\cdot3)\cdot 5/(5^3) = (16+9)\cdot5/125 = 125/125 = 1$. ✓

**Sigma 3 about [111] (60°):**

The Sigma 3 CSL is generated by a 60° rotation about the $[111]$ axis.
The exact rational matrix is:

$$R_{\Sigma 3} = \frac{1}{3}\begin{pmatrix} 2 & -1 & 2 \\ 2 & 2 & -1 \\ -1 & 2 & 2 \end{pmatrix}$$

Verification (exact, using `Fraction`):
- $\text{Tr}(R) = (2+2+2)/3 = 2$. For $\theta=60°$: $1 + 2\cos60° = 1+1 = 2$. ✓
- Axis: the fixed vector satisfies $R\hat{k} = \hat{k}$. For $\hat{k} = (1,1,1)/\sqrt{3}$:
  $R(1,1,1)^T = \frac{1}{3}(3,3,3)^T = (1,1,1)^T$. ✓
- $R^T R = I$ (exact). ✓
- $\det R = +1$ (exact). ✓

> **Previous error:** The matrix $\frac{1}{3}[[1,-2,2],[2,1,-2],[-2,2,1]]$ has
> $\text{Tr} = 1$ (implying $\theta = 90°$) and axis $[0,1,1]$, not $[111]$.
> The matrix $\frac{1}{3}[[-1,2,2],[2,-1,2],[2,2,-1]]$ is a 180° rotation about
> $[111]$ (also generating $\Sigma 3$, but as a twin relation, not a 60° twist).
> The correct 60° rotation matrix is the one given above.

---

## §12 Sigma (Coincidence Site Lattice)

### §12.1 Definition

For a rotation $R$ of the crystal lattice, the **Coincidence Site Lattice** (CSL)
is the set of lattice sites common to both the original and rotated lattices:

$$\text{CSL} = L \cap R(L)$$

**Sigma** ($\Sigma$) is defined as the ratio:

$$\Sigma = \frac{V_\text{CSL cell}}{V_\text{crystal cell}}$$

For a cubic crystal, $\Sigma$ is the smallest positive integer such that $\Sigma R$
is an integer matrix. For cubic symmetry, we additionally take the **odd part**
(divide out all factors of 2) since cubic symmetry relates $\Sigma$ to $\Sigma/2$ for
the same grain boundary type.

### §12.2 Algorithm

Given rotation matrix $R$ (with float or exact rational entries):

1. Convert each entry to a `Fraction` (if float, use `Fraction(entry).limit_denominator(10**9)`).
2. Find $q = \text{lcm}$ of all denominators of the nine `Fraction` entries.
3. Form the integer matrix $A = qR$ (each entry is now an exact integer).
4. Find $g = \gcd$ of all entries of $A$ (using `math.gcd` iteratively).
5. Set $\Sigma_0 = q / g$ (this is a positive integer by construction).
6. **Remove the 2-adic factor:** divide out the largest power of 2 dividing $\Sigma_0$:

$$\Sigma = \Sigma_0 \text{ with all factors of 2 divided out}$$

In code:
```python
while Sigma_0 % 2 == 0:
    Sigma_0 //= 2
Sigma = Sigma_0
```

> **Previous error:** Step 5 used $\Sigma_0 / 2^{\lfloor \log_2 \Sigma_0 \rfloor}$,
> which removes at most one factor of 2 (and yields non-integers like $6/4 = 1.5$
> when $\Sigma_0 = 6$). The correct operation is to divide out **all** factors of 2
> via the 2-adic valuation $\nu_2(\Sigma_0)$: $\Sigma = \Sigma_0 / 2^{\nu_2(\Sigma_0)}$.
> The while-loop above computes this exactly.

For exact CSL rotations (§11.1), this yields $\Sigma = 5$ and $\Sigma = 3$ exactly.

---

## §13 Bicrystal Construction

### §13.1 Geometry

- **Grain 1**: The reference crystal occupying the half-space $d \le 0$
  (signed distance $\le 0$ from the boundary plane).
- **Grain 2**: The rotated crystal occupying $d > 0$.
- **Boundary plane**: Defined in Grain 1's frame by $(hkl)$ via §3.
- **Rotation centre**: A point on the boundary plane.

### §13.2 Atom Positions

1. Generate a supercell large enough to include $N$ layers on each side.
2. For Grain 1: keep atoms with $d \le 0$.
3. For Grain 2: apply $R$ to all atom positions, then keep atoms with $d > 0$.
4. Remove from Grain 2 any atom within distance $\delta_\text{thresh}$ of any
   Grain 1 atom, where:
   $$\delta_\text{thresh} = \alpha \cdot d_\text{NN}$$
   with $\alpha \in (0, 1)$ (default $0.5$), and $d_\text{NN}$ from §10.

### §13.3 Boundary Classification

Let $\hat{k}$ be the rotation axis (unit vector) and $\hat{n}$ the boundary normal.

| Condition | Type |
|-----------|------|
| $|\hat{k} \cdot \hat{n}| < \epsilon$ | Tilt boundary |
| $||\hat{k} \cdot \hat{n}| - 1| < \epsilon$ | Twist boundary |
| Otherwise | Mixed boundary |

For exact CSL presets use exact arithmetic for this test.

### §13.4 Projection Modes

Two distinct projection modes are supported, both using the shared projection
engine from §6:

#### §13.4a Physical Bicrystal Mode (Half-Space)

The standard bicrystal construction from §13.1–13.2:
- Grain 1 (reference): atoms with $d \le 0$.
- Grain 2 (rotated): atoms with $d > 0$, after removal of overlapping atoms.
- Plan view: both grains projected onto the boundary plane via §6.2.
- The CSL 2D cell area is $\Sigma \times A_{2D}$ of the reference grain's 2D cell.
- The fraction of coincident sites in the plan view is $1/\Sigma$.

#### §13.4b Dichromatic Pattern Mode (Interpenetrating)

In the **dichromatic pattern** (also called the **dichromatic complex**), both
Grain 1 and the rotated Grain 2 span the **full** spherical/cylindrical aperture
(§1.4) **without** half-space truncation and **without** atom removal:
- All atoms of Grain 1 within the aperture are shown.
- All atoms of the rotated Grain 2 within the aperture are shown (including the
  region $d \le 0$).
- Coincident sites — positions common to both lattices within a tolerance
  $\varepsilon_\text{CSL}$ (default: $0.05 a$) — are identified and highlighted.

This mode displays the **full CSL pattern** and is used to visualise coincident
site fractions and the DSC lattice.

### §13.5 Visual Encoding Convention

The marker encoding is standardised to match crystallographic projection slides:

| Mode | Attribute mapped to shape | Attribute mapped to fill/colour |
|------|--------------------------|----------------------------------|
| **Single-crystal projection** | — (all circles) | Layer index (A/B/C/… → distinct colours) |
| **Bicrystal / dichromatic** | Stacking label → shape (see palette below) | Grain identity (Grain 1 → filled black, Grain 2 → open/blue, Coincident CSL → highlighted red) |

**Shape palette** (matplotlib marker codes, cycled modulo the palette length for
high-index planes with more than 4 distinct stacking layers):

```python
SHAPE_PALETTE = ['o', 's', '^', 'D', 'p', 'h', 'v', '<', '>']
# circle, square, triangle_up, diamond, pentagon, hexagon,
# triangle_down, triangle_left, triangle_right
```

Stacking label index $k$ (0-based: A=0, B=1, C=2, …) maps to
`SHAPE_PALETTE[k % len(SHAPE_PALETTE)]`. This ensures planes such as $(210)$
or $(311)$ with 6 or more distinct layers render without palette exhaustion.

For single-crystal projections, unique $(x_{2D}, y_{2D})$ sites are additionally
annotated with their fractional $(x, y, z)$ coordinates (after deduplication
from §1.4).

### §13.6 Side View (Tilt Boundaries)

Project onto the plane perpendicular to the tilt axis. The boundary appears as
a vertical line. Both grains are visible on either side. Atoms are coloured by
grain identity using the convention in §13.5.

---

## §14 Property-Based Test Invariants

The following must hold for **arbitrary** valid $(hkl)$:

1. $\mathbf{e}_1 \cdot \mathbf{e}_2 = 0$, $|\mathbf{e}_1| = |\mathbf{e}_2| = 1$,
   $\mathbf{e}_1 \cdot \hat{n} = \mathbf{e}_2 \cdot \hat{n} = 0$,
   $\mathbf{e}_1 \times \mathbf{e}_2 = \hat{n}$.

2. All **3D projected positions** lie in the $(hkl)$ plane:
   $$\hat{n} \cdot \mathbf{r}_\perp^{3D} = 0$$
   where $\mathbf{r}_\perp^{3D} = (\mathbf{r}\cdot\mathbf{e}_1)\mathbf{e}_1
   + (\mathbf{r}\cdot\mathbf{e}_2)\mathbf{e}_2$ is the 3D in-plane component of $\mathbf{r}$.
   (Note: $\mathbf{r}_\perp$ as a **2D** vector has no $\hat{n}$ component by
   construction; the invariant is on the 3D form $\mathbf{r}_\perp^{3D}$, which
   is orthogonal to $\hat{n}$ if and only if the projection is correct.)

3. Layer assignments are **exact** integers; no two atoms in the same layer
   disagree on their signed-distance key (computed via §4 exact arithmetic).

4. $\rho_{hkl} \times A_{2D} = 1$ exactly (a positive integer equal to 1) for
   SC, FCC, and BCC, since the 2D primitive cell derived from the true Bravais
   in-plane vectors (§7.1) contains exactly one atom. The test asserts
   **exact integer equality**, not approximate.

5. $(hkl)$ and $(-h,-k,-l)$ give identical $\rho_{hkl}$ and $d_{hkl}^\text{true}$,
   and the same **stacking period length** $P$ (the number of distinct layers per
   repeat period). The raw fractional-coordinate pairs $(s_1^\text{mod}, s_2^\text{mod})$
   are **not** required to be the same set: negating $\hat{n}$ reverses the layer
   ordering and may change which layer is indexed as 0, transforming the
   $(s_1, s_2)$ pairs by a unimodular shift and/or reordering.

   The **canonical** assertion is:
   - Let $\text{seq}(hkl)$ be the canonically relabelled stacking string of $(hkl)$
     (first-appearance labelling over one full period, sorted by ascending $d$).
   - Let $\text{seq}(-h,-k,-l)$ be the same for $(-h,-k,-l)$ (sorted by ascending
     $d$ under $-\hat{n}$, equivalently descending $d$ under $\hat{n}$).
   - Then $\text{seq}(-h,-k,-l)$ must equal the **canonical relabelling of the
     reversal** of $\text{seq}(hkl)$.

   Example: FCC (111) gives $\text{seq} = \text{ABC}$; reversed = $\text{CBA}$;
   relabelled (C→A, B→B, A→C) = $\text{ABC}$. So both orientations give
   $\text{ABC}$. The test asserts $P_{(hkl)} = P_{(-h,-k,-l)}$ and that the
   canonical-reversed string of $(hkl)$ equals the canonical string of $(-h,-k,-l)$.

6. Any $(hkl)$ related by the cubic point group (permutations and sign flips of
   indices, e.g., $(100)$, $(010)$, $(001)$, $(-1\,0\,0)$, etc.) gives the same
   $\rho_{hkl}$ and $d_{hkl}^\text{true}$ to within $10^{-12}$.
