# B-Matrix Method Comparison

This diagnostic harness compares the existing `finite_difference_B` baseline with the additive hybrid `analytical_B` API. It does not switch Stage 3D, Wilson GF, or VEDA-like diagnostics to the hybrid matrix.

Default four-file probe:

    .\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --out outputs\bmatrix_compare_minimal

Full fixture sweep:

    .\.venv312\Scripts\python.exe benchmarks\bmatrix_compare\compare_bmatrix_methods.py --full-sweep --out outputs\bmatrix_compare_full

Outputs:

- `bmatrix_method_comparison_summary.json`
- `bmatrix_method_comparison_summary.csv`
- `bmatrix_method_comparison_rows.csv`
- `bmatrix_method_comparison_selected_basis_differences.csv`

The summary reports per-file row deltas, hybrid row method counts, fallback reasons, redundant B rank/condition, independently selected basis rank/condition for each method, and whether selected basis indices changed at the same rank. Near-linear and high-angle rows are expected to fall back through `singular_or_near_linear_angle`.

The row CSV includes atom indices plus angle geometry for three-atom coordinates. The selected-basis-differences CSV records the finite-difference and hybrid rows chosen at each differing selected-basis position, replacement rank/condition/min-singular diagnostics, and whether replacing each swapped row preserves the original selected rank.

Current GAP 2 acceptance policy treats selected-basis index differences as allowed only when row-level deltas are within tolerance, redundant rank is unchanged, selected rank is unchanged, and every selected-basis replacement preserves rank. The current full sweep has 8 rank-preserving selected-basis differences in aromatic fixtures; these are diagnostic alternatives, not a production switch.

A `PASS` claim requires an actual command transcript; missing or unrun comparisons are not evidence.
