"""
Color and marker palettes for crystal plane and bicrystal visualization.

Implements visual encoding standards per MATH_SPEC §13.5 and FEATURES.md.
"""
from __future__ import annotations

# Marker shape palette per MATH_SPEC §13.5
# Cycled modulo palette length for high-index planes
SHAPE_PALETTE: list[str] = ["o", "s", "^", "D", "p", "h", "v", "<", ">"]

# Distinct layer color palette (accessible, high-contrast)
LAYER_COLOR_PALETTE: list[str] = [
    "#1f77b4",  # A - blue
    "#ff7f0e",  # B - orange
    "#2ca02c",  # C - green
    "#d62728",  # D - red
    "#9467bd",  # E - purple
    "#8c564b",  # F - brown
    "#e377c2",  # G - pink
    "#7f7f7f",  # H - gray
    "#bcbd22",  # I - yellow-green
    "#17becf",  # J - cyan
]

# Bicrystal grain visual encoding (MATH_SPEC §13.5)
GRAIN1_COLOR: str = "#111111"       # Grain 1: filled black/dark
GRAIN2_COLOR: str = "#1f77b4"       # Grain 2: open blue
COINCIDENT_COLOR: str = "#d62728"   # Coincident CSL: highlighted red


def get_shape_for_label(label: str) -> str:
    """Return marker shape for stacking label ('A'->'o', 'B'->'s', etc.)."""
    idx = ord(label.upper()[0]) - ord("A") if label else 0
    return SHAPE_PALETTE[idx % len(SHAPE_PALETTE)]


def get_color_for_layer(layer_idx: int) -> str:
    """Return distinct hex color for integer layer index without exhaustion."""
    return LAYER_COLOR_PALETTE[abs(layer_idx) % len(LAYER_COLOR_PALETTE)]
