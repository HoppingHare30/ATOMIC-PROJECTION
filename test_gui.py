"""
Headless GUI verification tests for Lattice Generator (Phase 4).

TP-73 through TP-77.

The tests run only when PySide6 is importable and can create a QApplication
without hanging. On macOS with the system Python 3.9 (no offscreen plugin
compiled), or any CI environment that lacks a display and has no offscreen
plugin, the entire module is skipped cleanly.
"""
from __future__ import annotations

import os
import sys
import pytest

# ---------------------------------------------------------------------------
# 1. Check PySide6 availability before importing anything from it
# ---------------------------------------------------------------------------
pyside6_available = True
pyside6_skip_reason = ""

try:
    # Import guard only – don't create a QApplication yet
    import importlib.util
    if importlib.util.find_spec("PySide6") is None:
        pyside6_available = False
        pyside6_skip_reason = "PySide6 not installed"
except Exception as e:
    pyside6_available = False
    pyside6_skip_reason = f"PySide6 import error: {e}"

# ---------------------------------------------------------------------------
# 2. Try a quick canary: can we actually create a QApplication offscreen?
# ---------------------------------------------------------------------------
if pyside6_available:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        from PySide6.QtWidgets import QApplication
        from PySide6.QtCore import QTimer

        _canary_app = QApplication.instance()
        if _canary_app is None:
            _canary_app = QApplication.__new__(QApplication)
            # Attempt creation with a 2-second abort so the test suite never hangs
            import signal

            def _timeout(signum, frame):
                raise RuntimeError("QApplication creation timed out")

            old_handler = signal.signal(signal.SIGALRM, _timeout)
            signal.alarm(5)
            try:
                _canary_app = QApplication(sys.argv[:1])
            except Exception as e:
                pyside6_available = False
                pyside6_skip_reason = f"QApplication creation failed: {e}"
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
    except Exception as e:
        pyside6_available = False
        pyside6_skip_reason = f"PySide6 Qt init failed: {e}"

skip_if_no_pyside = pytest.mark.skipif(
    not pyside6_available,
    reason=pyside6_skip_reason or "PySide6 unavailable or display missing",
)

# ---------------------------------------------------------------------------
# 3. Only import the MainWindow if PySide6 is confirmed working
# ---------------------------------------------------------------------------
if pyside6_available:
    from lattice_generator.gui.app import MainWindow


# ---------------------------------------------------------------------------
# 4. Session-scoped QApplication fixture
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication instance (skipped if PySide6 unavailable)."""
    if not pyside6_available:
        pytest.skip(pyside6_skip_reason)
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv[:1])
    yield app


# ---------------------------------------------------------------------------
# 5. Tests (TP-73 – TP-77)
# ---------------------------------------------------------------------------
@skip_if_no_pyside
def test_gui_single_crystal_mode(qapp):
    """TP-73 — MainWindow initialisation and FCC (111) single-crystal metrics."""
    window = MainWindow()
    assert window.windowTitle().startswith("Lattice Generator")

    window.structure_combo.setCurrentText("fcc")
    window.h_spin.setValue(1)
    window.k_spin.setValue(1)
    window.l_spin.setValue(1)
    window.n_radius_spin.setValue(2)

    window.recompute()

    metrics = window.metrics_text.toPlainText()
    assert "FCC" in metrics
    assert "(1 1 1)" in metrics
    assert "Planar Density" in metrics
    assert "Stacking Sequence" in metrics

    res = window._current_result
    assert res["mode"] == "single_crystal"
    assert res["structure"] == "fcc"
    assert res["hkl"] == [1, 1, 1]
    assert len(res["positions"]) > 0

    window.close()


@skip_if_no_pyside
def test_gui_bicrystal_csl_mode(qapp):
    """TP-74 / TP-75 — Bicrystal mode with CSL Sigma 5 preset."""
    window = MainWindow()

    window.mode_combo.setCurrentIndex(1)
    assert window.bicrystal_box.isVisible()

    window.structure_combo.setCurrentText("fcc")
    window.h_spin.setValue(0)
    window.k_spin.setValue(0)
    window.l_spin.setValue(1)
    window.n_radius_spin.setValue(2)

    window.rot_mode_combo.setCurrentIndex(0)
    window.csl_preset_combo.setCurrentText("sigma5_001")

    window.recompute()

    metrics = window.metrics_text.toPlainText()
    assert "CSL Index (Sigma): \u03a35" in metrics
    assert "Misorientation Angle" in metrics

    res = window._current_result
    assert res["mode"] == "bicrystal"
    assert res["sigma"] == 5
    assert res["removed_count"] >= 0

    window.close()


@skip_if_no_pyside
def test_gui_bicrystal_custom_rotation(qapp):
    """TP-76 — Bicrystal with custom rotation axis and angle."""
    window = MainWindow()

    window.mode_combo.setCurrentIndex(1)
    window.rot_mode_combo.setCurrentIndex(1)  # Custom mode

    window.axis_x_spin.setValue(1.0)
    window.axis_y_spin.setValue(0.0)
    window.axis_z_spin.setValue(0.0)
    window.angle_spin.setValue(90.0)

    window.recompute()

    res = window._current_result
    assert res["mode"] == "bicrystal"
    assert abs(res["misorientation_deg"] - 90.0) < 0.1

    window.close()
