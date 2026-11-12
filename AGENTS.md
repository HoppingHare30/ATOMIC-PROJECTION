# AGENTS.md — Agent Rules and Conventions

> This document codifies the rules that govern all agent (AI-assisted) work
> on this project. It is the operational constitution: stable across all phases.
> Source-of-truth for scope: SOURCE_OF_TRUTH.md.

---

## §1 Non-Negotiable Priorities

1. **Mathematical correctness above all.** If speed, looks, or convenience
   conflict with correctness, correctness wins. No exceptions.
2. **Never loosen test tolerances to make a test pass.** Find and fix the
   mathematical cause. Document the fix in the phase's status log
   (SOURCE_OF_TRUTH.md §4).
3. **Never edit expected test values to match wrong output** without first
   proving from MATH_SPEC.md that the expected value was wrong.

---

## §2 Phase Gate Protocol

Work proceeds in strictly ordered phases. At the end of each phase:

1. Run the full test suite: `pytest tests/ -v`.
2. Report the results (number passed, failed, skipped) in this conversation.
3. If any gate test fails, fix it before declaring the phase complete.
4. Write a short summary in SOURCE_OF_TRUTH.md §4 (Status Log) with:
   - Phase number and date.
   - Tests run and results.
   - Any issues found and their root causes.
5. **Stop and wait for explicit user approval** before beginning the next phase.

---

## §3 Code Conventions

### §3.1 Docstrings

Every `core/` function **must** have a docstring with:

```python
def function_name(args) -> return_type:
    """
    One-line description.

    Implements: <equation or reference to MATH_SPEC §N>.
    Units: <input and output units, e.g. 'positions in units of a'>
    Exact: <whether exact arithmetic is used and why>

    Args:
        arg_name: description
    Returns:
        description
    Raises:
        ExceptionType: when
    """
```

### §3.2 Arithmetic Policy (from SOURCE_OF_TRUTH.md §2.2)

| Quantity | Required type | Rationale |
|----------|--------------|-----------|
| Layer index / plane membership | `fractions.Fraction` or doubled-int | MATH_SPEC §4 |
| Stacking offsets | `fractions.Fraction` | MATH_SPEC §9.2 |
| Sigma, CSL | Exact integer arithmetic | MATH_SPEC §12 |
| CSL rotation matrix entries | `fractions.Fraction` | MATH_SPEC §11.1 |
| Arbitrary rotation matrix | `float` (numpy) | MATH_SPEC §11 |
| Rendered positions, bond lengths | `float` | Cast only at viz boundary |

### §3.3 Import Boundaries

```
core/     → may import: numpy, fractions, math, itertools, typing
           must NOT import: matplotlib, PySide6, any viz/ or gui/ module
viz/      → may import: matplotlib, numpy, any core/ module
           must NOT import: PySide6, any gui/ module
export/   → may import: numpy, json, csv, any core/ module
           must NOT import: matplotlib (except for PNG/SVG save), PySide6
gui/      → may import: PySide6, any core/, viz/, export/ module
           must NOT import: numpy (except to read arrays for display), fractions
tests/    → may import: pytest, hypothesis, numpy, any core/ module
           must NOT import: matplotlib, PySide6 (except in GUI tests)
```

### §3.4 Naming Conventions

| Entity | Convention | Example |
|--------|-----------|---------|
| Module | lowercase_underscore | `cell2d.py` |
| Public function | lowercase_underscore | `assign_layers()` |
| Private helper | _prefixed | `_gauss_step()` |
| Constant | UPPER_CASE | `FCC_BASIS` |
| Test function | `test_<description>` | `test_fcc_111_density()` |
| Test hypothesis strategy | `@given(...)` | `@given(hkl_strategy())` |

### §3.5 Units

- All internal positions: **units of `a`** (lattice parameter = 1).
- All energies, forces: not applicable (this is a geometry tool).
- All areas: **units of `a²`** (hence planar density in `a⁻²`).
- All angles: **degrees** in user-facing output; **radians** in numpy trig calls.

---

## §4 Testing Rules

1. Tests live in `tests/` and import from `core/` only.
2. Every reference value in TEST_PLAN.md is implemented as an automated test.
3. Hypothesis property-based tests run for ≥ 100 examples (default setting).
4. Tolerance is `1e-12` unless a looser tolerance is explicitly justified
   with a comment citing the mathematical reason.
5. Integer quantities (layer index, CN, Sigma) must pass with **exact equality**,
   never approximate.
6. After any code change, the **full** test suite is run, not just the
   affected tests.

---

## §5 Bicrystal-Specific Rules

1. The bicrystal builder calls `core/projection.py` for all geometry; it does
   not duplicate projection code.
2. CSL preset rotation matrices are stored as `fractions.Fraction` arrays; they
   are never computed from floating-point trigonometry.
3. `R^T R = I` and `det R = +1` are **asserted** in `build_bicrystal()`, not
   just tested.
4. The atom-removal threshold is derived from `d_NN` returned by
   `core/neighbors.nearest_neighbour()`. It is never a hardcoded number.
5. Sigma is computed by the algorithm in MATH_SPEC §12.2 and must yield exact
   integers for all CSL presets.

---

## §6 Documentation Rules

1. SOURCE_OF_TRUTH.md is updated at the **start** of each phase (status log
   entry) and whenever a convention, decision, or reference value changes.
2. MATH_SPEC.md is updated if a derivation is found to be incomplete or wrong.
   Any change must be accompanied by a note explaining what was wrong and why.
3. FEATURES.md acceptance criteria are never weakened; if an AC is found to be
   wrong, it is updated with a note and the change is discussed with the user.
4. TEST_PLAN.md reference values are never changed without a proof from
   MATH_SPEC.md that the old value was wrong.

---

## §7 Communication Protocol

After completing each phase the agent will produce:

```
## Phase N Complete

### Tests run
<pytest summary output>

### Results
- Passed: N
- Failed: 0
- Skipped: 0

### Notes
<Any issues found, root causes, design decisions made>

### Status log entry (for SOURCE_OF_TRUTH.md)
<Draft of the §4 entry to be added>

Awaiting approval to proceed to Phase N+1.
```

If any test fails, the communication is:

```
## Phase N — BLOCKED (test failure)

### Failing tests
<pytest -v output for failing tests>

### Root cause analysis
<Mathematical diagnosis from MATH_SPEC.md>

### Fix applied
<What was changed and why>

### Re-run results
<pytest summary after fix>
```

---

## §8 Forbidden Patterns

The following patterns are forbidden and their presence is a defect:

| Pattern | Why forbidden |
|---------|--------------|
| `if abs(d1 - d2) < 1e-6: same_layer = True` | Float comparison for layer grouping |
| `expected = 2.3094` (hardcoded in test) | Use `4/math.sqrt(3)` or `import math; math.isclose(...)` |
| `tolerance = 1e-6` loosened to `1e-4` | Never loosen; fix the code |
| `import matplotlib` in `core/` | Breaks headless testing |
| `import PySide6` in `core/` or `viz/` | Breaks server-side use |
| `d_thresh = 0.3` (hardcoded) | Must derive from `d_NN` |
| Geometric computations in `gui/app.py` | GUI is a thin shell only |
| CSL rotation from `np.sin`, `np.cos` | Use exact `Fraction` entries |
