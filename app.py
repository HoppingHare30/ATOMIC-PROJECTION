"""
Main graphical user interface for Lattice Generator.

Thin GUI shell integrating core geometry calculation engines, visualization
plotting routines, and data export handlers.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys
from typing import Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from lattice_generator.core.bicrystal import (
    build_bicrystal,
    coincidence_cell_area,
    coincident_site_fraction,
    get_csl_preset,
    rodrigues_rotation,
)
from lattice_generator.core.cell2d import cell_area, find_inplane_vectors
from lattice_generator.core.density import (
    formula_spacing as interplanar_spacing_formula,
    planar_density,
    true_spacing as true_interplanar_spacing,
)
from lattice_generator.core.lattice import generate_supercell
from lattice_generator.core.layers import (
    assign_layers,
    filter_by_n_layers as filter_layers_by_n,
    filter_by_distance as filter_layers_by_distance,
    stacking_sequence,
)
from lattice_generator.core.neighbors import nearest_neighbour
from lattice_generator.core.plane import unit_normal, validate_and_reduce
from lattice_generator.core.projection import orthonormal_basis, project_2d
from lattice_generator.export.writer import (
    export_csv,
    export_image,
    export_json,
    export_xyz,
)
from lattice_generator.viz.plot_2d import plot_plan_view, plot_side_view
from lattice_generator.viz.plot_3d import plot_3d_structure


class MainWindow(QMainWindow):
    """Main window for crystal plane and bicrystal analysis."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Lattice Generator — Crystal Plane & Bicrystal Analyzer")
        self.resize(1280, 850)

        # State storage for exports
        self._current_result: dict[str, Any] = {}

        self._init_ui()
        self._update_visibility()
        self.recompute()

    def _init_ui(self) -> None:
        main_splitter = QSplitter(Qt.Horizontal)
        self.setCentralWidget(main_splitter)

        # ------------------- Left Controls Panel -------------------
        control_scroll = QScrollArea()
        control_scroll.setWidgetResizable(True)
        control_scroll.setMinimumWidth(380)
        control_scroll.setMaximumWidth(450)

        controls_container = QWidget()
        controls_layout = QVBoxLayout(controls_container)
        controls_layout.setSpacing(12)

        # Mode selector
        mode_box = QGroupBox("Analysis Mode")
        mode_layout = QVBoxLayout(mode_box)
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Single Crystal Plane", "Bicrystal Interface"])
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        mode_layout.addWidget(self.mode_combo)
        controls_layout.addWidget(mode_box)

        # Primary crystal parameters
        param_box = QGroupBox("Lattice & Plane Parameters")
        param_form = QFormLayout(param_box)

        self.structure_combo = QComboBox()
        self.structure_combo.addItems(["fcc", "bcc", "sc"])
        param_form.addRow("Structure:", self.structure_combo)

        # Miller indices
        hkl_layout = QHBoxLayout()
        self.h_spin = QSpinBox()
        self.h_spin.setRange(-20, 20)
        self.h_spin.setValue(1)
        self.k_spin = QSpinBox()
        self.k_spin.setRange(-20, 20)
        self.k_spin.setValue(1)
        self.l_spin = QSpinBox()
        self.l_spin.setRange(-20, 20)
        self.l_spin.setValue(1)
        hkl_layout.addWidget(QLabel("h:"))
        hkl_layout.addWidget(self.h_spin)
        hkl_layout.addWidget(QLabel("k:"))
        hkl_layout.addWidget(self.k_spin)
        hkl_layout.addWidget(QLabel("l:"))
        hkl_layout.addWidget(self.l_spin)
        param_form.addRow("Miller (h k l):", hkl_layout)

        # Supercell radius
        self.n_radius_spin = QSpinBox()
        self.n_radius_spin.setRange(1, 8)
        self.n_radius_spin.setValue(3)
        param_form.addRow("Supercell N:", self.n_radius_spin)

        # Layer filtering
        self.filter_mode_combo = QComboBox()
        self.filter_mode_combo.addItems(["All Layers", "Closest N Layers", "Distance Cutoff"])
        self.filter_mode_combo.currentIndexChanged.connect(self._on_filter_mode_changed)
        param_form.addRow("Layer Filter:", self.filter_mode_combo)

        self.filter_val_spin = QDoubleSpinBox()
        self.filter_val_spin.setRange(0.1, 50.0)
        self.filter_val_spin.setValue(3.0)
        self.filter_val_spin.setEnabled(False)
        self.filter_val_label = QLabel("Filter Val:")
        param_form.addRow(self.filter_val_label, self.filter_val_spin)

        controls_layout.addWidget(param_box)

        # Bicrystal parameters
        self.bicrystal_box = QGroupBox("Bicrystal Misorientation")
        b_form = QFormLayout(self.bicrystal_box)

        self.rot_mode_combo = QComboBox()
        self.rot_mode_combo.addItems(["CSL Preset", "Custom Axis & Angle"])
        self.rot_mode_combo.currentIndexChanged.connect(self._on_rot_mode_changed)
        b_form.addRow("Rotation Mode:", self.rot_mode_combo)

        self.csl_preset_combo = QComboBox()
        self.csl_preset_combo.addItems([
            "sigma3_111",
            "sigma5_001",
            "sigma5_100",
            "sigma7_111",
            "sigma9_110",
            "sigma11_110",
        ])
        b_form.addRow("CSL Preset:", self.csl_preset_combo)

        # Custom axis
        axis_layout = QHBoxLayout()
        self.axis_x_spin = QDoubleSpinBox()
        self.axis_x_spin.setRange(-10.0, 10.0)
        self.axis_x_spin.setValue(0.0)
        self.axis_y_spin = QDoubleSpinBox()
        self.axis_y_spin.setRange(-10.0, 10.0)
        self.axis_y_spin.setValue(0.0)
        self.axis_z_spin = QDoubleSpinBox()
        self.axis_z_spin.setRange(-10.0, 10.0)
        self.axis_z_spin.setValue(1.0)
        axis_layout.addWidget(self.axis_x_spin)
        axis_layout.addWidget(self.axis_y_spin)
        axis_layout.addWidget(self.axis_z_spin)
        b_form.addRow("Axis [u v w]:", axis_layout)

        self.angle_spin = QDoubleSpinBox()
        self.angle_spin.setRange(-360.0, 360.0)
        self.angle_spin.setValue(36.87)
        b_form.addRow("Angle (deg):", self.angle_spin)

        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.1, 1.0)
        self.alpha_spin.setSingleStep(0.05)
        self.alpha_spin.setValue(0.5)
        b_form.addRow("Overlap alpha (d_NN):", self.alpha_spin)

        controls_layout.addWidget(self.bicrystal_box)

        # Compute button
        self.compute_btn = QPushButton("Compute & Update View")
        self.compute_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.compute_btn.clicked.connect(self.recompute)
        controls_layout.addWidget(self.compute_btn)

        # Metrics display
        metrics_box = QGroupBox("Calculated Metrics")
        metrics_layout = QVBoxLayout(metrics_box)
        self.metrics_text = QTextEdit()
        self.metrics_text.setReadOnly(True)
        self.metrics_text.setMinimumHeight(180)
        metrics_layout.addWidget(self.metrics_text)
        controls_layout.addWidget(metrics_box)

        # Export group
        export_box = QGroupBox("Export Data & Plots")
        export_layout = QVBoxLayout(export_box)

        btn_row1 = QHBoxLayout()
        self.export_csv_btn = QPushButton("CSV")
        self.export_csv_btn.clicked.connect(self._export_csv)
        self.export_xyz_btn = QPushButton("XYZ")
        self.export_xyz_btn.clicked.connect(self._export_xyz)
        btn_row1.addWidget(self.export_csv_btn)
        btn_row1.addWidget(self.export_xyz_btn)
        export_layout.addLayout(btn_row1)

        btn_row2 = QHBoxLayout()
        self.export_json_btn = QPushButton("JSON")
        self.export_json_btn.clicked.connect(self._export_json)
        self.export_plot_btn = QPushButton("Plot (PNG/SVG)")
        self.export_plot_btn.clicked.connect(self._export_plot)
        btn_row2.addWidget(self.export_json_btn)
        btn_row2.addWidget(self.export_plot_btn)
        export_layout.addLayout(btn_row2)

        controls_layout.addWidget(export_box)
        controls_layout.addStretch()

        control_scroll.setWidget(controls_container)
        main_splitter.addWidget(control_scroll)

        # ------------------- Right Visualization Panel -------------------
        self.tabs = QTabWidget()

        # Tab 1: Plan View
        self.tab_plan = QWidget()
        plan_layout = QVBoxLayout(self.tab_plan)
        self.fig_plan = Figure(figsize=(7, 6))
        self.canvas_plan = FigureCanvasQTAgg(self.fig_plan)
        self.toolbar_plan = NavigationToolbar2QT(self.canvas_plan, self.tab_plan)
        plan_layout.addWidget(self.toolbar_plan)
        plan_layout.addWidget(self.canvas_plan)
        self.tabs.addTab(self.tab_plan, "2D Plan View")

        # Tab 2: Side View
        self.tab_side = QWidget()
        side_layout = QVBoxLayout(self.tab_side)
        self.fig_side = Figure(figsize=(7, 5))
        self.canvas_side = FigureCanvasQTAgg(self.fig_side)
        self.toolbar_side = NavigationToolbar2QT(self.canvas_side, self.tab_side)
        side_layout.addWidget(self.toolbar_side)
        side_layout.addWidget(self.canvas_side)
        self.tabs.addTab(self.tab_side, "Side View (Tilt View)")

        # Tab 3: 3D Structure
        self.tab_3d = QWidget()
        layout_3d = QVBoxLayout(self.tab_3d)
        self.fig_3d = Figure(figsize=(7, 6))
        self.canvas_3d = FigureCanvasQTAgg(self.fig_3d)
        self.toolbar_3d = NavigationToolbar2QT(self.canvas_3d, self.tab_3d)
        layout_3d.addWidget(self.toolbar_3d)
        layout_3d.addWidget(self.canvas_3d)
        self.tabs.addTab(self.tab_3d, "3D Structure")

        main_splitter.addWidget(self.tabs)
        main_splitter.setStretchFactor(1, 2)

    def _on_mode_changed(self) -> None:
        self._update_visibility()
        self.recompute()

    def _on_filter_mode_changed(self) -> None:
        mode_idx = self.filter_mode_combo.currentIndex()
        if mode_idx == 0:
            self.filter_val_spin.setEnabled(False)
        elif mode_idx == 1:
            self.filter_val_spin.setEnabled(True)
            self.filter_val_label.setText("Layer count:")
            self.filter_val_spin.setDecimals(0)
            self.filter_val_spin.setValue(4.0)
        else:
            self.filter_val_spin.setEnabled(True)
            self.filter_val_label.setText("Distance cutoff:")
            self.filter_val_spin.setDecimals(2)
            self.filter_val_spin.setValue(2.0)

    def _on_rot_mode_changed(self) -> None:
        is_preset = self.rot_mode_combo.currentIndex() == 0
        self.csl_preset_combo.setEnabled(is_preset)
        self.axis_x_spin.setEnabled(not is_preset)
        self.axis_y_spin.setEnabled(not is_preset)
        self.axis_z_spin.setEnabled(not is_preset)
        self.angle_spin.setEnabled(not is_preset)

    def _update_visibility(self) -> None:
        is_bicrystal = self.mode_combo.currentIndex() == 1
        self.bicrystal_box.setVisible(is_bicrystal)
        self._on_rot_mode_changed()

    def recompute(self) -> None:
        """Trigger calculation and refresh all visualizations and metrics."""
        h, k, l = self.h_spin.value(), self.k_spin.value(), self.l_spin.value()
        if h == 0 and k == 0 and l == 0:
            QMessageBox.critical(self, "Invalid Input", "Miller indices (h, k, l) cannot all be zero.")
            return

        structure = self.structure_combo.currentText().lower()
        n_radius = self.n_radius_spin.value()
        is_bicrystal = self.mode_combo.currentIndex() == 1

        try:
            h_red, k_red, l_red = validate_and_reduce(h, k, l)
            n_hat = unit_normal(h_red, k_red, l_red)
            e1, e2, _ = orthonormal_basis(h_red, k_red, l_red)

            # Nearest neighbour distance and coordination number
            d_nn, cn = nearest_neighbour(structure)
            # Spacing and density
            d_form = interplanar_spacing_formula(h_red, k_red, l_red)
            d_true = true_interplanar_spacing(structure, h_red, k_red, l_red)
            rho = planar_density(structure, h_red, k_red, l_red)
            v1_3d, v2_3d = find_inplane_vectors(structure, h_red, k_red, l_red)
            v1_2d = project_2d(v1_3d, e1, e2)
            v2_2d = project_2d(v2_3d, e1, e2)
            a_2d = cell_area(v1_3d, v2_3d)

            if not is_bicrystal:
                # ----------------- Single Crystal Mode -----------------
                atoms = generate_supercell(structure, n_radius)
                layer_dict = assign_layers(atoms, h_red, k_red, l_red)

                # Filter layers if requested
                filt_idx = self.filter_mode_combo.currentIndex()
                if filt_idx == 1:
                    layer_dict = filter_layers_by_n(
                        layer_dict,
                        atoms,
                        h_red, k_red, l_red,
                        int(self.filter_val_spin.value()),
                    )
                elif filt_idx == 2:
                    layer_dict = filter_layers_by_distance(
                        layer_dict,
                        atoms,
                        h_red, k_red, l_red,
                        0.0,
                        self.filter_val_spin.value(),
                    )

                # Collect atoms that survived filtering
                selected_atom_indices: list[int] = []
                for indices in layer_dict.values():
                    selected_atom_indices.extend(indices)
                selected_atom_indices.sort()

                if not selected_atom_indices:
                    QMessageBox.warning(self, "Empty Slab", "No atoms found in the selected layer filter range.")
                    return

                filtered_atoms = atoms[selected_atom_indices]
                # Reassign layer keys for filtered subset
                filtered_layers = assign_layers(filtered_atoms, h_red, k_red, l_red)
                layer_keys_flat = [0] * len(filtered_atoms)
                for k_val, idxs in filtered_layers.items():
                    for idx in idxs:
                        layer_keys_flat[idx] = k_val

                stacking_labels = stacking_sequence(
                    filtered_layers,
                    filtered_atoms,
                    v1_3d,
                    v2_3d,
                )

                # Project coordinates
                proj_coords = project_2d(filtered_atoms, e1, e2)
                u = proj_coords[:, 0]
                v = proj_coords[:, 1]
                # Perpendicular signed distances
                d_signed = [float(sum(filtered_atoms[i, j] * n_hat[j] for j in range(3))) for i in range(len(filtered_atoms))]

                # Clear and render figures
                self.fig_plan.clear()
                ax_plan = self.fig_plan.add_subplot(111)
                plot_plan_view(
                    u,
                    v,
                    layer_keys_flat,
                    stacking_labels=stacking_labels,
                    v1_2d=v1_2d,
                    v2_2d=v2_2d,
                    title=f"Plan View — {structure.upper()} ({h_red} {k_red} {l_red})",
                    ax=ax_plan,
                )
                self.canvas_plan.draw()

                self.fig_side.clear()
                ax_side = self.fig_side.add_subplot(111)
                plot_side_view(
                    u,
                    d_signed,
                    title=f"Side View — {structure.upper()} ({h_red} {k_red} {l_red})",
                    ax=ax_side,
                )
                self.canvas_side.draw()

                self.fig_3d.clear()
                ax_3d = self.fig_3d.add_subplot(111, projection="3d")
                plot_3d_structure(
                    filtered_atoms,
                    layer_keys=layer_keys_flat,
                    plane_normal=n_hat,
                    title=f"3D Structure — {structure.upper()} ({h_red} {k_red} {l_red})",
                    ax=ax_3d,
                )
                self.canvas_3d.draw()

                # Update metrics
                metrics_lines = [
                    f"<b>Structure:</b> {structure.upper()}",
                    f"<b>Miller Plane:</b> ({h_red} {k_red} {l_red})",
                    f"<b>Formula Spacing d_hkl:</b> {d_form:.6f} a",
                    f"<b>True Interplanar Spacing:</b> {d_true:.6f} a",
                    f"<b>Planar Density:</b> {rho:.6f} a⁻²",
                    f"<b>2D Unit Cell Area:</b> {a_2d:.6f} a²",
                    f"<b>2D Basis v1:</b> [{v1_3d[0]}, {v1_3d[1]}, {v1_3d[2]}]",
                    f"<b>2D Basis v2:</b> [{v2_3d[0]}, {v2_3d[1]}, {v2_3d[2]}]",
                    f"<b>Stacking Sequence:</b> {''.join(stacking_labels[:min(12, len(stacking_labels))])}...",
                    f"<b>Nearest Neighbour d_NN:</b> {d_nn:.6f} a",
                    f"<b>Coordination Number (CN):</b> {cn}",
                    f"<b>Total Atoms in Slab:</b> {len(filtered_atoms)}",
                ]
                self.metrics_text.setHtml("<br>".join(metrics_lines))

                # Store result
                self._current_result = {
                    "mode": "single_crystal",
                    "structure": structure,
                    "hkl": [h_red, k_red, l_red],
                    "positions": filtered_atoms,
                    "layer_keys": layer_keys_flat,
                    "stacking_labels": stacking_labels,
                    "planar_density": rho,
                    "interplanar_spacing": d_true,
                }

            else:
                # ----------------- Bicrystal Mode -----------------
                rot_mode = self.rot_mode_combo.currentIndex()
                if rot_mode == 0:
                    preset_name = self.csl_preset_combo.currentText()
                    R = get_csl_preset(preset_name)
                else:
                    axis_vec = (
                        self.axis_x_spin.value(),
                        self.axis_y_spin.value(),
                        self.axis_z_spin.value(),
                    )
                    theta_deg = self.angle_spin.value()
                    R = rodrigues_rotation(axis_vec, math.radians(theta_deg))

                alpha_val = self.alpha_spin.value()
                b_res = build_bicrystal(
                    structure,
                    h_red,
                    k_red,
                    l_red,
                    R,
                    N=n_radius,
                    removal_alpha=alpha_val,
                )

                positions = b_res["positions"]
                grain_ids = b_res["grain_ids"]
                layer_keys = b_res["layer_keys"]

                proj_coords = project_2d(positions, e1, e2)
                u = proj_coords[:, 0]
                v = proj_coords[:, 1]
                d_signed = [float(sum(positions[i, j] * n_hat[j] for j in range(3))) for i in range(len(positions))]

                # CSL cell area and site fraction
                sigma_val = b_res["sigma"]
                csl_area = coincidence_cell_area(sigma_val, a_2d)
                csl_fraction = coincident_site_fraction(sigma_val)

                # Clear and render figures
                self.fig_plan.clear()
                ax_plan = self.fig_plan.add_subplot(111)
                plot_plan_view(
                    u,
                    v,
                    layer_keys,
                    grain_ids=grain_ids,
                    v1_2d=v1_2d,
                    v2_2d=v2_2d,
                    title=f"Bicrystal Plan View — {structure.upper()} ({h_red} {k_red} {l_red}) Σ{sigma_val}",
                    ax=ax_plan,
                )
                self.canvas_plan.draw()

                self.fig_side.clear()
                ax_side = self.fig_side.add_subplot(111)
                plot_side_view(
                    u,
                    d_signed,
                    grain_ids=grain_ids,
                    title=f"Bicrystal Side View — Boundary {b_res['boundary_type'].upper()}",
                    ax=ax_side,
                )
                self.canvas_side.draw()

                self.fig_3d.clear()
                ax_3d = self.fig_3d.add_subplot(111, projection="3d")
                plot_3d_structure(
                    positions,
                    grain_ids=grain_ids,
                    plane_normal=n_hat,
                    title=f"3D Bicrystal — Σ{sigma_val} ({b_res['boundary_type']})",
                    ax=ax_3d,
                )
                self.canvas_3d.draw()

                # Update metrics
                metrics_lines = [
                    f"<b>Bicrystal Boundary:</b> {b_res['boundary_type'].upper()}",
                    f"<b>CSL Index (Sigma):</b> Σ{sigma_val}",
                    f"<b>Misorientation Angle:</b> {b_res['misorientation_deg']:.3f}°",
                    f"<b>Coincident Site Fraction (1/Σ):</b> {csl_fraction:.4f}",
                    f"<b>Coincidence Cell Area:</b> {csl_area:.6f} a²",
                    f"<b>Reference 2D Cell Area:</b> {a_2d:.6f} a²",
                    f"<b>Grain 1 Atoms (d ≤ 0):</b> {len(b_res['grain1_positions'])}",
                    f"<b>Grain 2 Atoms (d > 0):</b> {len(b_res['grain2_positions'])}",
                    f"<b>Removed Overlap Atoms:</b> {b_res['removed_count']}",
                    f"<b>Removal Threshold:</b> {alpha_val * d_nn:.6f} a ({alpha_val:.2f} × d_NN)",
                ]
                self.metrics_text.setHtml("<br>".join(metrics_lines))

                # Store result
                self._current_result = {
                    "mode": "bicrystal",
                    "structure": structure,
                    "hkl": [h_red, k_red, l_red],
                    "positions": positions,
                    "grain_ids": grain_ids,
                    "layer_keys": layer_keys,
                    "sigma": sigma_val,
                    "misorientation_deg": b_res["misorientation_deg"],
                    "boundary_type": b_res["boundary_type"],
                    "removed_count": b_res["removed_count"],
                }

        except Exception as e:
            QMessageBox.critical(self, "Computation Error", str(e))

    def _export_csv(self) -> None:
        if not self._current_result:
            return
        filename, _ = QFileDialog.getSaveFileName(self, "Export CSV", "", "CSV Files (*.csv)")
        if not filename:
            return
        extra_cols = {}
        if "grain_ids" in self._current_result:
            extra_cols["grain_id"] = self._current_result["grain_ids"]
        export_csv(
            filename,
            self._current_result["positions"],
            self._current_result["layer_keys"],
            labels=self._current_result.get("stacking_labels"),
            extra_cols=extra_cols if extra_cols else None,
        )
        QMessageBox.information(self, "Export Complete", f"Data exported successfully to:\n{filename}")

    def _export_xyz(self) -> None:
        if not self._current_result:
            return
        filename, _ = QFileDialog.getSaveFileName(self, "Export XYZ", "", "XYZ Files (*.xyz)")
        if not filename:
            return
        structure = self._current_result.get("structure", "fcc")
        species = "Cu" if structure == "fcc" else ("Fe" if structure == "bcc" else "Po")
        comment = f"{self._current_result.get('mode')} {structure} hkl={self._current_result.get('hkl')}"
        export_xyz(filename, self._current_result["positions"], species=species, comment=comment)
        QMessageBox.information(self, "Export Complete", f"XYZ file exported successfully to:\n{filename}")

    def _export_json(self) -> None:
        if not self._current_result:
            return
        filename, _ = QFileDialog.getSaveFileName(self, "Export JSON", "", "JSON Files (*.json)")
        if not filename:
            return
        export_json(filename, self._current_result)
        QMessageBox.information(self, "Export Complete", f"JSON metadata exported successfully to:\n{filename}")

    def _export_plot(self) -> None:
        filename, _ = QFileDialog.getSaveFileName(self, "Save Plot", "", "PNG Image (*.png);;SVG Vector (*.svg)")
        if not filename:
            return
        curr_tab = self.tabs.currentIndex()
        fig = self.fig_plan if curr_tab == 0 else (self.fig_side if curr_tab == 1 else self.fig_3d)
        export_image(fig, filename)
        QMessageBox.information(self, "Plot Saved", f"Plot exported successfully to:\n{filename}")


def launch_app() -> None:
    """Launch the Lattice Generator GUI application."""
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
