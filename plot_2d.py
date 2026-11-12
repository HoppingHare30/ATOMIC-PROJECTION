"""
2D plan view and side view plotting for crystal planes and bicrystals.

Implements visual encoding standards from MATH_SPEC §13.5, FEATURES.md F04/F12.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless testability
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

from lattice_generator.viz.palettes import (
    COINCIDENT_COLOR,
    GRAIN1_COLOR,
    GRAIN2_COLOR,
    get_color_for_layer,
    get_shape_for_label,
)


def plot_plan_view(
    u: np.ndarray,
    v: np.ndarray,
    layer_keys: list[int] | np.ndarray,
    stacking_labels: list[str] | None = None,
    grain_ids: np.ndarray | None = None,
    coincident_mask: np.ndarray | None = None,
    v1_2d: np.ndarray | None = None,
    v2_2d: np.ndarray | None = None,
    csl_v1_2d: np.ndarray | None = None,
    csl_v2_2d: np.ndarray | None = None,
    title: str = "In-Plane Projection (Plan View)",
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (8.0, 7.0),
) -> plt.Figure:
    """
    Generate 2D plan view of crystal planes or bicrystal interfaces.

    Args:
        u: 1D array of horizontal in-plane coordinates (units of a).
        v: 1D array of vertical in-plane coordinates (units of a).
        layer_keys: 1D array/list of layer indices.
        stacking_labels: Optional list of stacking labels ('A', 'B', 'C', ...).
        grain_ids: Optional 1D array of grain IDs (1 or 2). If None, treated as single crystal.
        coincident_mask: Optional boolean array marking coincident CSL sites.
        v1_2d: Optional 2D primitive vector 1 (u1, v1).
        v2_2d: Optional 2D primitive vector 2 (u2, v2).
        csl_v1_2d: Optional 2D CSL cell vector 1.
        csl_v2_2d: Optional 2D CSL cell vector 2.
        title: Plot title.
        ax: Optional Matplotlib Axes to draw into.
        figsize: Figure dimensions (width, height) in inches.

    Returns:
        Matplotlib Figure object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    u_arr = np.asarray(u, dtype=np.float64)
    v_arr = np.asarray(v, dtype=np.float64)
    n_atoms = len(u_arr)

    if grain_ids is None:
        # Single crystal mode: color by layer / stacking label
        unique_keys = sorted(list(set(layer_keys)))
        for i, key in enumerate(unique_keys):
            mask = np.array([k == key for k in layer_keys])
            label_text = f"Layer {key}"
            marker = "o"
            if stacking_labels is not None and len(stacking_labels) == len(unique_keys):
                st_label = stacking_labels[i]
                label_text = f"Layer {key} ({st_label})"
                marker = get_shape_for_label(st_label)

            color = get_color_for_layer(key)
            ax.scatter(
                u_arr[mask],
                v_arr[mask],
                c=color,
                marker=marker,
                s=60,
                alpha=0.85,
                edgecolors="black",
                linewidths=0.5,
                label=label_text,
                zorder=3,
            )
    else:
        # Bicrystal mode per MATH_SPEC §13.5
        g_ids = np.asarray(grain_ids, dtype=np.int64)

        # Grain 1: filled black markers
        m1 = g_ids == 1
        if np.any(m1):
            ax.scatter(
                u_arr[m1],
                v_arr[m1],
                c=GRAIN1_COLOR,
                marker="o",
                s=50,
                alpha=0.8,
                label="Grain 1 (Reference)",
                zorder=3,
            )

        # Grain 2: open blue markers
        m2 = g_ids == 2
        if np.any(m2):
            ax.scatter(
                u_arr[m2],
                v_arr[m2],
                facecolors="none",
                edgecolors=GRAIN2_COLOR,
                marker="o",
                s=60,
                linewidths=1.5,
                alpha=0.9,
                label="Grain 2 (Rotated)",
                zorder=4,
            )

        # Coincident CSL sites (if present)
        if coincident_mask is not None and np.any(coincident_mask):
            ax.scatter(
                u_arr[coincident_mask],
                v_arr[coincident_mask],
                c=COINCIDENT_COLOR,
                marker="s",
                s=90,
                edgecolors="black",
                linewidths=1.0,
                label="Coincident CSL Site",
                zorder=5,
            )

    # Draw 2D primitive unit cell parallelogram if provided
    if v1_2d is not None and v2_2d is not None:
        p0 = np.array([0.0, 0.0])
        p1 = np.asarray(v1_2d, dtype=np.float64)[:2]
        p2 = p1 + np.asarray(v2_2d, dtype=np.float64)[:2]
        p3 = np.asarray(v2_2d, dtype=np.float64)[:2]
        cell_poly = Polygon(
            [p0, p1, p2, p3],
            closed=True,
            fill=False,
            edgecolor="#2ca02c",
            linestyle="--",
            linewidth=2.0,
            label="2D Primitive Cell",
            zorder=2,
        )
        ax.add_patch(cell_poly)

    # Draw CSL 2D unit cell parallelogram if provided
    if csl_v1_2d is not None and csl_v2_2d is not None:
        c0 = np.array([0.0, 0.0])
        c1 = np.asarray(csl_v1_2d, dtype=np.float64)[:2]
        c2 = c1 + np.asarray(csl_v2_2d, dtype=np.float64)[:2]
        c3 = np.asarray(csl_v2_2d, dtype=np.float64)[:2]
        csl_poly = Polygon(
            [c0, c1, c2, c3],
            closed=True,
            fill=False,
            edgecolor="#d62728",
            linestyle="-",
            linewidth=2.0,
            label="CSL 2D Cell",
            zorder=2,
        )
        ax.add_patch(csl_poly)

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(r"In-Plane Coordinate $u$ ($a$)", fontsize=11)
    ax.set_ylabel(r"In-Plane Coordinate $v$ ($a$)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.5)

    if n_atoms > 0:
        ax.legend(loc="upper right", framealpha=0.9, fontsize=9)

    fig.tight_layout()
    return fig


def plot_side_view(
    u: np.ndarray,
    signed_distances: np.ndarray,
    grain_ids: np.ndarray | None = None,
    title: str = "Side View (Tilt View)",
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (8.0, 5.0),
) -> plt.Figure:
    """
    Generate side view (edge-on view) showing atomic layers and interface boundary.

    Args:
        u: 1D array of horizontal in-plane coordinates (units of a).
        signed_distances: 1D array of signed perpendicular distances d (units of a).
        grain_ids: Optional 1D array of grain IDs (1 or 2).
        title: Plot title.
        ax: Optional Matplotlib Axes.
        figsize: Figure dimensions in inches.

    Returns:
        Matplotlib Figure object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    u_arr = np.asarray(u, dtype=np.float64)
    d_arr = np.asarray(signed_distances, dtype=np.float64)

    if grain_ids is not None:
        g_ids = np.asarray(grain_ids, dtype=np.int64)
        m1 = g_ids == 1
        m2 = g_ids == 2
        if np.any(m1):
            ax.scatter(u_arr[m1], d_arr[m1], c=GRAIN1_COLOR, s=40, alpha=0.8, label=r"Grain 1 ($d \leq 0$)")
        if np.any(m2):
            ax.scatter(u_arr[m2], d_arr[m2], facecolors="none", edgecolors=GRAIN2_COLOR, s=40, linewidths=1.2, label=r"Grain 2 ($d > 0$)")
    else:
        ax.scatter(u_arr, d_arr, c="#1f77b4", s=40, alpha=0.8, edgecolors="black", linewidths=0.5)

    # Draw boundary plane line at d = 0
    if len(u_arr) > 0:
        min_u, max_u = float(np.min(u_arr)), float(np.max(u_arr))
        margin = max(0.5, (max_u - min_u) * 0.05)
        ax.plot(
            [min_u - margin, max_u + margin],
            [0.0, 0.0],
            color="#d62728",
            linestyle="--",
            linewidth=1.8,
            label="Boundary Plane ($d=0$)",
            zorder=2,
        )

    ax.set_xlabel(r"In-Plane Coordinate $u$ ($a$)", fontsize=11)
    ax.set_ylabel(r"Perpendicular Distance $d$ ($a$)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)

    fig.tight_layout()
    return fig
