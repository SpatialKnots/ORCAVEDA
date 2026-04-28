# ORCAVEDA Stage 3D v5.0 Regression Infrastructure Report

## Status

Regression harness status: **PASS**

Rows: 94  
PASS: 94  
FAIL: 0

## Scope

Stage 3D v5.0 formalizes regression testing. Chemistry and assignment behavior are inherited from v4.9.

## Regression Set

- Acetone_freq.hess
- CH3CN_freq.hess
- DMF_freq.hess
- DMSO_freq.hess
- EtOH_freq.hess
- MeOH_freq.hess
- NMP_freq.hess
- iPrOH_freq.hess

## Main Checks

- Expected molecule outputs present
- `rank_B_independent == expected_rank_3N_minus_6`
- Normal-mode orientation rule reported
- Expected functional groups detected
- Mid-frequency diagnostic assignments present
- Forbidden CH2 labels absent where applicable
- No high-frequency unassigned modes
- Monoethanolamine sanity checks skipped for external molecules
- Protected X-H unused norm/power cleanup preserved
- DMF negative-frequency warning tolerated as expected input-quality warning

## Remaining Warning

DMF retains an expected input-quality warning: `negative_freq_count_after_first_6 = 2`.

## Method Boundary

This remains a PED-like / assignment-audit layer, not a strict VEDA-equivalent PED implementation.
