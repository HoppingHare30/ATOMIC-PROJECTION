"""
3D crystal and bicrystal structure visualization using Matplotlib mplot3d.

Implements 3D atomic scatter plots and crystallographic plane highlights.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import numpy as np

from lattice_generator.viz.palettes import GRAIN1_COLOR, GRAIN2_COLOR, get_color_for_layer


def plot_3d_structure(
    positions: np.ndarray,
    grain_ids: np.ndarray | None = None,
    layer_keys: list[int] | np.ndarray | None = None,
    plane_normal: np.ndarray | None = None,
    plane_d: float = 0.0,
    title: str = "3D Crystal Structure",
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (8.0, 7.0),
) -> plt.Figure:
    """
    Render 3D scatter plot of atomic positions with optional boundary plane highlight.

    Args:
        positions: Nx3 array of atom coordinates in units of a.
        grain_ids: Optional 1D array of grain IDs (1 or 2).
        layer_keys: Optional 1D array of discrete layer indices.
        plane_normal: Optional 3D plane unit normal for drawing a plane surface.
        plane_d: Distance offset for plane surface.
        title: Plot title.
        ax: Optional 3D Matplotlib Axes.
        figsize: Figure dimensions in inches.

    Returns:
        Matplotlib Figure object.
    """
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.figure

    pos = np.asarray(positions, dtype=np.float64)
    if len(pos) == 0:
        return fig

    xs, ys, zs = pos[:, 0], pos[:, 1], pos[:, 2]

    if grain_ids is not None:
        g_ids = np.asarray(grain_ids, dtype=np.int64)
        m1 = g_ids == 1
        m2 = g_ids == 2
        if np.any(m1):
            ax.scatter(xs[m1], ys[m1], zs[m1], c=GRAIN1_COLOR, s=35, alpha=0.85, label="Grain 1")
        if np.any(m2):
            ax.scatter(xs[m2], ys[m2], zs[m2], c=GRAIN2_COLOR, s=35, alpha=0.85, label="Grain 2")
    elif layer_keys is not None:
        unique_keys = sorted(list(set(layer_keys)))
        for key in unique_keys:
            mask = np.array([k == key for k in layer_keys])
            color = get_color_for_layer(key)
            ax.scatter(xs[mask], ys[mask], zs[mask], c=color, s=35, alpha=0.85, label=f"Layer {key}")
    else:
        ax.scatter(xs, ys, zs, c="#1f77b4", s=35, alpha=0.85, edgecolors="black", linewidths=0.3)

    # Optional plane rendering
    if plane_normal is not None:
        n = np.asarray(plane_normal, dtype=np.float64) / np.linalg.norm(plane_normal)
        min_x, max_x = float(np.min(xs)), float(np.max(xs))
        min_y, max_y = float(np.min(ys)), float(np.max(ys))
        grid_x, grid_y = np.meshgrid(
            np.linspace(min_x, max_x, 10),
            np.linspace(min_y, max_y, 10),
        )
        if abs(n[2]) > 1e-6:
            grid_z = (plane_d - n[0] * grid_x - n[1] * grid_y) / n[2]
            ax.plot_surface(grid_x, grid_y, grid_z, alpha=0.25, color="red")

    ax.set_xlabel(r"$x / a$", fontsize=10)
    ax.set_ylabel(r"$y / a$", fontsize=10)
    ax.set_zlabel(r"$z / a$", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    return fig
