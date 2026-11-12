"""
Bicrystal generation, Rodrigues rotations, CSL presets, and Sigma computation.

Constructs bicrystal interface slabs from two oriented crystal grains,
computes Coincidence Site Lattice (CSL) Sigma factors using exact rational arithmetic,
and handles boundary classification and atom overlap removal.
"""
from __future__ import annotations

from fractions import Fraction
import math
import numpy as np

from lattice_generator.core.lattice import generate_supercell
from lattice_generator.core.plane import unit_normal, validate_and_reduce
from lattice_generator.core.neighbors import nearest_neighbour


def rodrigues_rotation(
    axis: tuple[float, float, float] | np.ndarray,
    theta_rad: float,
) -> np.ndarray:
    """
    Construct a 3x3 rotation matrix using Rodrigues' rotation formula.

    Implements: MATH_SPEC §11:
      K = skew(k_hat)
      R = I*cos(theta) + (1 - cos(theta))*(k_hat outer k_hat) + sin(theta)*K
    Units: Axis dimensionless direction vector; theta_rad in radians.
    Exact: Computed in float64 with rigorous assertion R^T R = I and det R = +1.

    Args:
        axis: Rotation axis as 3-tuple or NumPy array of length 3.
        theta_rad: Rotation angle in radians.

    Returns:
        3x3 NumPy array (float64) representing the proper rotation matrix.

    Raises:
        ValueError: If axis is zero vector.
        AssertionError: If R^T R != I or det(R) != +1 within 1e-12.
    """
    ax = np.asarray(axis, dtype=np.float64)
    norm = np.linalg.norm(ax)
    if norm < 1e-12:
        raise ValueError("Rotation axis cannot be zero vector.")
    k = ax / norm

    cos_t = math.cos(theta_rad)
    sin_t = math.sin(theta_rad)

    # Skew-symmetric cross-product matrix K
    k0, k1, k2 = k[0], k[1], k[2]
    K = np.array([
        [0.0, -k2, k1],
        [k2, 0.0, -k0],
        [-k1, k0, 0.0],
    ], dtype=np.float64)

    # Rodrigues formula
    I = np.eye(3, dtype=np.float64)
    R = I * cos_t + (1.0 - cos_t) * np.outer(k, k) + sin_t * K

    # Invariant assertions (AGENTS.md §5 rule 3)
    assert np.allclose(R.T @ R, I, atol=1e-12), "Rodrigues rotation is not orthogonal."
    assert abs(np.linalg.det(R) - 1.0) < 1e-12, "Rodrigues rotation determinant is not +1."

    return R


def get_csl_preset(name: str) -> np.ndarray:
    """
    Return exact rational 3x3 rotation matrix for standard CSL presets.

    Implements: MATH_SPEC §11.1, §12.3:
      sigma5_001: 36.87° about [001]
      sigma5_100: 36.87° about [100]
      sigma3_111: 60° about [111]
      sigma7_111: 38.21° / 81.79° about [111]
      sigma9_110: 38.94° / 141.06° about [110]
      sigma11_110: 50.48° about [110]
    Units: Dimensionless rotation matrix entries in Q.
    Exact: All matrix entries stored as fractions.Fraction objects.

    Args:
        name: Identifier string for preset (e.g. 'sigma3_111', 'sigma5_001').

    Returns:
        3x3 NumPy array of dtype=object containing Fraction instances.

    Raises:
        ValueError: If preset name is not recognized.
    """
    key = name.strip().lower()
    if key in ("sigma5_001", "sigma5_[001]", "sigma5_z"):
        # MATH_SPEC §11.1
        mat = [
            [Fraction(4, 5), Fraction(-3, 5), Fraction(0, 1)],
            [Fraction(3, 5), Fraction(4, 5), Fraction(0, 1)],
            [Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)],
        ]
    elif key in ("sigma5_100", "sigma5_[100]", "sigma5_x"):
        # MATH_SPEC §11.1
        mat = [
            [Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)],
            [Fraction(0, 1), Fraction(4, 5), Fraction(-3, 5)],
            [Fraction(0, 1), Fraction(3, 5), Fraction(4, 5)],
        ]
    elif key in ("sigma3_111", "sigma3_[111]"):
        # MATH_SPEC §11.1
        mat = [
            [Fraction(2, 3), Fraction(-1, 3), Fraction(2, 3)],
            [Fraction(2, 3), Fraction(2, 3), Fraction(-1, 3)],
            [Fraction(-1, 3), Fraction(2, 3), Fraction(2, 3)],
        ]
    elif key in ("sigma7_111", "sigma7_[111]"):
        mat = [
            [Fraction(3, 7), Fraction(-2, 7), Fraction(6, 7)],
            [Fraction(6, 7), Fraction(3, 7), Fraction(-2, 7)],
            [Fraction(-2, 7), Fraction(6, 7), Fraction(3, 7)],
        ]
    elif key in ("sigma9_110", "sigma9_[110]"):
        mat = [
            [Fraction(1, 9), Fraction(8, 9), Fraction(4, 9)],
            [Fraction(8, 9), Fraction(1, 9), Fraction(-4, 9)],
            [Fraction(-4, 9), Fraction(4, 9), Fraction(-7, 9)],
        ]
    elif key in ("sigma11_110", "sigma11_[110]"):
        mat = [
            [Fraction(9, 11), Fraction(2, 11), Fraction(6, 11)],
            [Fraction(2, 11), Fraction(9, 11), Fraction(-6, 11)],
            [Fraction(-6, 11), Fraction(6, 11), Fraction(7, 11)],
        ]
    else:
        raise ValueError(
            f"Unknown CSL preset '{name}'. Supported: 'sigma3_111', 'sigma5_001', "
            f"'sigma5_100', 'sigma7_111', 'sigma9_110', 'sigma11_110'."
        )

    return np.array(mat, dtype=object)


def compute_sigma(R: np.ndarray) -> int:
    """
    Compute exact integer Sigma (CSL index) from rotation matrix R.

    Implements: MATH_SPEC §12.2:
      1. Convert entries of R to Fraction.
      2. q = lcm of all denominators.
      3. A = q * R (integer matrix).
      4. g = gcd of all entries of A.
      5. Sigma_0 = q // g.
      6. Remove 2-adic factor: while Sigma_0 % 2 == 0: Sigma_0 //= 2.
    Units: Dimensionless integer ratio.
    Exact: Exact integer arithmetic without floating-point approximations.

    Args:
        R: 3x3 rotation matrix (object array of Fractions or float NumPy array).

    Returns:
        Exact integer Sigma value.
    """
    fractions_list: list[Fraction] = []
    for r in range(3):
        for c in range(3):
            val = R[r, c]
            if isinstance(val, Fraction):
                fractions_list.append(val)
            else:
                fractions_list.append(Fraction(float(val)).limit_denominator(10000))

    # LCM of denominators
    q = 1
    for f in fractions_list:
        q = math.lcm(q, f.denominator)

    # Integer matrix A
    a_entries = [int(f * q) for f in fractions_list]

    # GCD of all entries of A
    g = abs(a_entries[0])
    for entry in a_entries[1:]:
        g = math.gcd(g, abs(entry))

    sigma_0 = q // g if g != 0 else q

    # Remove all factors of 2 per MATH_SPEC §12.2
    while sigma_0 % 2 == 0:
        sigma_0 //= 2

    return int(sigma_0)


def extract_rotation_axis_angle(
    R: np.ndarray,
) -> tuple[np.ndarray, float]:
    """
    Extract unit rotation axis and angle from rotation matrix R.

    Implements: MATH_SPEC §11, §13.3.
    Units: Axis is unit vector; angle returned in radians.
    Exact: Float64 trigonometric inverse.

    Args:
        R: 3x3 rotation matrix.

    Returns:
        Tuple (axis, angle_rad) where axis has shape (3,) and |axis| = 1.0.
    """
    r_mat = np.asarray(R, dtype=np.float64)
    tr = np.trace(r_mat)
    cos_t = float(np.clip((tr - 1.0) / 2.0, -1.0, 1.0))
    theta = math.acos(cos_t)

    if theta < 1e-12:
        return (np.array([0.0, 0.0, 1.0], dtype=np.float64), 0.0)

    if abs(theta - math.pi) < 1e-12:
        # At theta = pi: R + I = 2 * outer(k, k)
        m = r_mat + np.eye(3, dtype=np.float64)
        diag = np.diag(m)
        best_col = int(np.argmax(diag))
        col = m[:, best_col]
        norm = np.linalg.norm(col)
        axis = col / norm if norm > 1e-12 else np.array([1.0, 0.0, 0.0], dtype=np.float64)
        return (axis, math.pi)

    # General case
    k = np.array([
        r_mat[2, 1] - r_mat[1, 2],
        r_mat[0, 2] - r_mat[2, 0],
        r_mat[1, 0] - r_mat[0, 1],
    ], dtype=np.float64)
    norm = np.linalg.norm(k)
    axis = k / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0], dtype=np.float64)
    return (axis, theta)


def classify_boundary(
    axis: np.ndarray,
    normal: np.ndarray,
    tol: float = 1e-5,
) -> str:
    """
    Classify grain boundary as 'tilt', 'twist', or 'mixed'.

    Implements: MATH_SPEC §13.3:
      axis perp normal (|axis.normal| < tol) -> 'tilt'
      axis || normal (||axis.normal| - 1| < tol) -> 'twist'
      otherwise -> 'mixed'
    Units: Dimensionless vectors.
    Exact: Dot product compared with tolerance.

    Args:
        axis: Rotation axis unit vector.
        normal: Boundary plane normal unit vector.
        tol: Geometric alignment tolerance (default: 1e-5).

    Returns:
        String: 'tilt', 'twist', or 'mixed'.
    """
    k = np.asarray(axis, dtype=np.float64) / np.linalg.norm(axis)
    n = np.asarray(normal, dtype=np.float64) / np.linalg.norm(normal)
    dot = abs(float(np.dot(k, n)))

    if dot < tol:
        return "tilt"
    elif abs(dot - 1.0) < tol:
        return "twist"
    else:
        return "mixed"


def coincidence_cell_area(sigma: int, a_2d: float) -> float:
    """
    Compute coincidence cell area in plan view: A_CSL = Sigma * A_2D.

    Implements: MATH_SPEC §13.4a, FEATURES F12-AC12.
    Units: Area in units of a^2.
    Exact: Product of exact Sigma and 2D unit cell area.

    Args:
        sigma: Exact integer CSL Sigma factor.
        a_2d: Area of 2D primitive unit cell in units of a^2.

    Returns:
        Float coincidence cell area.
    """
    return float(sigma) * a_2d


def coincident_site_fraction(sigma: int) -> float:
    """
    Compute fraction of coincident sites in plan view: 1 / Sigma.

    Implements: MATH_SPEC §13.4a, FEATURES F12-AC17.
    Units: Dimensionless fraction.
    Exact: Reciprocal of exact integer Sigma.

    Args:
        sigma: Exact integer CSL Sigma factor.

    Returns:
        Float value 1.0 / sigma.
    """
    if sigma < 1:
        raise ValueError(f"Sigma must be positive integer, got {sigma}")
    return 1.0 / float(sigma)


class BicrystalResult(dict):
    """
    Structured result of bicrystal construction.

    Acts as both a dictionary and a 3-tuple (positions, grain_ids, layer_keys)
    for seamless unpacking and flexible inspection.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.__dict__.update(kwargs)

    def __iter__(self):
        yield self["positions"]
        yield self["grain_ids"]
        yield self["layer_keys"]

    def __getitem__(self, key):
        if isinstance(key, int):
            return [self["positions"], self["grain_ids"], self["layer_keys"]][key]
        return super().__getitem__(key)


def build_bicrystal(
    structure: str,
    h: int,
    k: int,
    l: int,
    R: np.ndarray,
    N: int = 3,
    removal_alpha: float = 0.5,
) -> BicrystalResult:
    """
    Construct bicrystal interface slab with overlap removal.

    Implements: MATH_SPEC §13.1–13.3:
      - Grain 1: d <= 0 (unrotated reference).
      - Grain 2: d > 0 (rotated by R).
      - Atom removal: Grain 2 atoms within alpha * d_NN of any Grain 1 atom are removed.
    Units: Coordinates in units of lattice parameter a.
    Exact: Assert R^T R = I and det R = +1 on entry; overlap threshold derived from d_NN.

    Args:
        structure: Crystal structure identifier ('fcc', 'bcc', 'sc').
        h: Miller index h of boundary plane.
        k: Miller index k of boundary plane.
        l: Miller index l of boundary plane.
        R: 3x3 rotation matrix (NumPy array or Fraction array).
        N: Supercell radius for centered slab (default: 3).
        removal_alpha: Fraction of d_NN for overlap threshold (default: 0.5).

    Returns:
        BicrystalResult containing:
          - 'positions': ndarray shape (M, 3) of combined atomic positions.
          - 'grain_ids': ndarray shape (M,) of grain identifiers (1 or 2).
          - 'layer_keys': ndarray shape (M,) of discrete layer indices.
          - 'grain1_positions': ndarray shape (M1, 3).
          - 'grain2_positions': ndarray shape (M2, 3).
          - 'removed_count': int count of removed overlapping atoms.
          - 'boundary_type': str ('tilt', 'twist', or 'mixed').
          - 'sigma': int exact CSL Sigma index.
          - 'misorientation_deg': float misorientation angle in degrees.

    Raises:
        AssertionError: If R is not a valid proper rotation matrix.
    """
    # 1. Assert orthogonality and determinant +1 on entry (AGENTS.md §5 rule 3)
    r_float = np.asarray(R, dtype=np.float64)
    assert np.allclose(r_float.T @ r_float, np.eye(3), atol=1e-12), "R^T R != I on entry to build_bicrystal"
    assert abs(np.linalg.det(r_float) - 1.0) < 1e-12, "det(R) != +1 on entry to build_bicrystal"

    # 2. Reference geometry
    h_red, k_red, l_red = validate_and_reduce(h, k, l)
    n_hat = unit_normal(h_red, k_red, l_red)

    # 3. Supercell atoms in centered box
    raw_atoms = generate_supercell(structure, N)

    # Check for 0° rotation (single crystal mode per TP-59)
    axis, theta_rad = extract_rotation_axis_angle(r_float)
    misorientation_deg = math.degrees(theta_rad)
    is_identity = (theta_rad < 1e-10) or np.allclose(r_float, np.eye(3), atol=1e-10)

    if is_identity:
        # TP-59: all atoms in one grain, removed_count = 0
        grain1_pos = raw_atoms
        grain2_pos = np.empty((0, 3), dtype=np.float64)
        positions = raw_atoms
        grain_ids = np.ones(len(raw_atoms), dtype=np.int64)
        removed_count = 0
    else:
        # 4. Grain 1: unrotated atoms with d <= 0
        d1 = np.dot(raw_atoms, n_hat)
        grain1_mask = d1 <= 1e-12
        grain1_pos = raw_atoms[grain1_mask]

        # 5. Grain 2: rotated atoms with d > 0
        # r' = R @ r
        rotated_all = (r_float @ raw_atoms.T).T
        d2 = np.dot(rotated_all, n_hat)
        grain2_mask = d2 > 1e-12
        grain2_pos_initial = rotated_all[grain2_mask]

        # 6. Overlap removal threshold derived from d_NN (AGENTS.md §5 rule 4)
        d_nn, _ = nearest_neighbour(structure)
        delta_thresh = removal_alpha * d_nn

        # Check pairwise distance to Grain 1 atoms
        grain2_kept: list[np.ndarray] = []
        removed_count = 0

        if len(grain1_pos) > 0 and len(grain2_pos_initial) > 0:
            strict_thresh = delta_thresh - 1e-12
            for g2_atom in grain2_pos_initial:
                dists = np.linalg.norm(grain1_pos - g2_atom, axis=1)
                min_dist = float(np.min(dists))
                if min_dist < strict_thresh:
                    removed_count += 1
                else:
                    grain2_kept.append(g2_atom)
        else:
            grain2_kept = list(grain2_pos_initial)

        grain2_pos = np.array(grain2_kept, dtype=np.float64) if grain2_kept else np.empty((0, 3), dtype=np.float64)

        n1 = len(grain1_pos)
        n2 = len(grain2_pos)
        positions = np.vstack([grain1_pos, grain2_pos]) if n2 > 0 else grain1_pos
        grain_ids = np.array([1] * n1 + [2] * n2, dtype=np.int64)

    # 8. Layer keys
    # For Grain 1: exact integer key = 2*(hx+ky+lz)
    fh, fk, fl = float(h_red), float(k_red), float(l_red)
    layer_keys_list: list[int] = []
    for p in grain1_pos:
        val = 2.0 * (fh * p[0] + fk * p[1] + fl * p[2])
        layer_keys_list.append(int(round(val)))

    # For Grain 2: integer layer keys along n_hat
    norm_hkl = math.sqrt(h_red * h_red + k_red * k_red + l_red * l_red)
    for p in grain2_pos:
        dist_n = float(np.dot(p, n_hat))
        val = 2.0 * norm_hkl * dist_n
        layer_keys_list.append(int(round(val)))

    layer_keys = np.array(layer_keys_list, dtype=np.int64)

    # 9. Misorientation angle and boundary classification
    axis, theta_rad = extract_rotation_axis_angle(r_float)
    misorientation_deg = math.degrees(theta_rad)
    b_type = classify_boundary(axis, n_hat)
    sigma_val = compute_sigma(R)

    return BicrystalResult(
        positions=positions,
        grain_ids=grain_ids,
        layer_keys=layer_keys,
        grain1_positions=grain1_pos,
        grain2_positions=grain2_pos,
        removed_count=removed_count,
        boundary_type=b_type,
        sigma=sigma_val,
        misorientation_deg=misorientation_deg,
    )
