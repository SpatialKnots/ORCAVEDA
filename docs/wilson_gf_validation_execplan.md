# Wilson GF Diagonalization Validation Prototype

This ExecPlan is a living document. The sections `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` must be kept up to date as work proceeds.

This plan follows `PLANS.md` at the repository root.

## Purpose / Big Picture

ORCAVEDA currently provides Stage 3D, a geometric and weighted independent-coordinate assignment audit, plus PED diagnostic layers that expose internal-coordinate contributions. The existing Wilson-style PED audit uses a Wilson `G = B M^-1 B^T` matrix and a reconstructed internal force matrix, but it does not solve and validate a closed Wilson GF eigenvalue problem.

This plan adds an opt-in Wilson GF diagonalization validation prototype. After implementation, a user should be able to run ORCAVEDA with `--wilson-gf-validation` and inspect new CSV outputs showing whether a selected nonredundant internal-coordinate basis can reproduce ORCA vibrational frequencies through a GF eigenvalue calculation. The initial output is a validation and diagnostic artifact, not a production final-label engine and not a claim of VEDA equivalence.

The first scientific goal is modest and testable: for H2O, verify whether the prototype can build a full-rank `3N-6` internal-coordinate basis, solve the GF eigenproblem, map the resulting vibrational eigenvalues to ORCA positive vibrational modes, and report fixed-conversion and empirical-ratio diagnostics without unsupported claims. Only after H2O passes should the same path be evaluated on NH3, formaldehyde, and ethene.

## Progress

- [x] (2026-05-14) Initial ExecPlan written from the reviewed Wilson GF upgrade proposal.
- [x] (2026-05-14) Phase 0 H2O numerical gate recorded: `data\hess\H2O_freq.hess` has 3 atoms, `$hessian` present, 3 positive ORCA vibrational modes, 3 internal coordinates selected, B-rank 3, and B-condition 110.8498941006996.
- [x] (2026-05-14) Phase 1 standalone prototype added in `src\wilson_gf.py`, plus focused tests in `tests\test_wilson_gf.py`.
- [x] (2026-05-14) Phase 2 H2O fixed-conversion validation passes after replacing the Euclidean PED pseudoinverse force reconstruction with a mass-metric Wilson GF internal-force reconstruction in the standalone prototype. H2O max relative error is `1.818626413830463e-05`.
- [x] (2026-05-14) Phase 3 opt-in pipeline integration added behind `--wilson-gf-validation`. Default runs do not include Wilson GF validation tables; opt-in H2O run writes validation, closed-PED, and basis-diagnostics CSVs.
- [x] (2026-05-14) Phase 4 target validation runs executed for H2O, NH3, formaldehyde, and ethene. An initial ethene run exposed a weak default basis (`G-rank 11`, `G-condition 2.5804969975545938e12`); after adding a small-system Wilson-GF-conditioned basis fallback, all four target molecules report `PASS`.
- [x] (2026-05-14) Broader opt-in validation batch run for MeOH, EtOH, acetone, CH3CN, DMSO, acetamide, phenol, and benzene. An initial CH3CN run reported `WARN`; after adding opt-in Wilson GF linear-bend components for near-linear bends, CH3CN reports `PASS` with `linear_bend_coordinate_used`.
- [x] (2026-05-14) Larger-molecule opt-in validation expansion run for pyridine, aniline, benzonitrile, nitrobenzene, benzoic acid, and acetophenone. All six CLI runs exited 0 and report `PASS`; aniline carries explicit positive-mode-count warnings.
- [x] (2026-05-14) Heavy opt-in validation batch run for NMP, N-methylaniline, piperidine, cyclohexane chair, acetanilide, and three monoethanolamine dimers. An initial acetanilide run reported `FAIL` due to G-rank loss, severe ill-conditioning, and positive-mode count below expected rank.
- [x] (2026-05-14) Added a deterministic mass-weighted residual-pivot validation fallback for large systems where exhaustive conditioned-basis search is infeasible. Acetanilide, acetophenone, and nitrobenzene now use better-conditioned validation bases and report `PASS`.
- [x] (2026-05-14) Extended large-system fallback acceptance to consider F-condition as well as G-condition. N-methylaniline now uses a better-conditioned validation basis and no longer carries `f_ill_conditioned`.
- [x] (2026-05-14) Hardened positive-mode-count diagnostics. Aniline, N-methylaniline, and acetanilide each have six zero ORCA modes plus one negative ORCA mode, so the opt-in Wilson GF CSVs now report nonpositive ORCA/GF counts and minimum nonpositive values.

## Surprises & Discoveries

- Observation: Current ORCAVEDA already has a Wilson GF-style PED audit in `src/ped.py`.
  Evidence: `build_wilson_g_matrix`, `reconstruct_internal_force_matrix`, and `build_wilson_ped_audit_dataframe` exist in `src/ped.py`; the audit projects ORCA normal modes through the selected internal-coordinate basis and reports G/F diagnostics, but it does not solve a closed GF eigenvalue problem.

- Observation: A nonredundant internal-coordinate GF matrix with `m = 3N-6` or `m = 3N-5` cannot produce six translational/rotational zero eigenvalues.
  Evidence: `G`, `F`, and `GF` in such a basis are `m x m`, so the eigenproblem yields `m` vibrational eigenvalues. Translational/rotational modes are removed by the nonredundant internal-coordinate basis and must not be expected as six zero columns in this prototype.

- Observation: The current parser and assignment code depend on the ORCA normal-mode orientation `normal_modes[:, mode]`.
  Evidence: `PLANS.md`, `src/ped.py`, and `src/ORCAVEDA_patched_stage3D_v5_0.py` explicitly use and protect this orientation.

- Observation: H2O fixed SI conversion failed when the first prototype reused the existing PED audit Euclidean pseudoinverse force reconstruction, but passed when the standalone GF prototype used the mass-metric Wilson back-transform `J = M^-1 B^T G^-1`.
  Evidence: Euclidean pseudoinverse reconstruction on `data\hess\H2O_freq.hess` gave reconstructed frequencies `1540.560643, 4049.375290, 4131.042557` cm-1 and max relative error `0.002300580590651947`. Mass-metric reconstruction gave reconstructed frequencies `1540.708842, 4049.448994, 4140.643569` cm-1 and max relative error `1.818626413830463e-05`.

- Observation: Ethene required a Wilson-GF-conditioned basis fallback; the default independent-coordinate selector was full rank in B but too ill-conditioned for closed GF validation.
  Evidence: Before the fallback, ethene selected indices `0;1;2;3;4;9;10;11;12;13;14;15` with B-rank `12` and B-condition about `3.146e6`, but Wilson G-rank `11`, G-condition `2580496997554.5938`, F-condition `5156767183685.332`, and warnings `basis_rank_below_expected; g_ill_conditioned; f_ill_conditioned`. The fallback chooses `0;4;5;6;7;9;12;13;14;15;17;18`, giving G-rank `12`, G-condition `28.556444875327514`, F-rank `12`, F-condition `65.64477191436387`, and no warnings.

- Observation: CH3CN is the first broader-batch molecule where fixed conversion fails despite full G/F rank.
  Evidence: `outputs\wilson_gf_batch_ch3cn\CH3CN__wilson_gf_validation.csv` reports `WARN`, basis size `12`, expected rank `12`, G-rank `12`, G-condition `423673515.5779151`, F-rank `12`, F-condition `4979156243.108764`, max relative error `0.018454913619246567`, and warnings `near_linear_bend_coordinate; fixed_conversion_failed; empirical_ratio_only`. The largest errors are in mode 6 (`389.6377536209096` cm-1, reconstructed `382.44702253503846` cm-1, relative error `0.018454913619246567`) and mode 7 (`390.57271369641705` cm-1, reconstructed `389.46928984530246` cm-1, relative error `0.0028251432125703077`).

- Observation: Exhaustive primitive-basis search did not fix CH3CN.
  Evidence: An exhaustive diagnostic over all 1136 full-rank 12-coordinate primitive subsets for `data\hess\CH3CN_freq.hess` found zero subsets with `PASS` status and `max_relative_error <= 1.0e-4`. The best observed max relative error was `0.01845477680255841`. The CH3CN geometry contains a near-linear `ang(C1-C2-N3)` coordinate at `179.9940447723222` degrees, and the low-frequency warning modes are dominated by this near-linear bend/torsion subspace.

- Observation: Opt-in Wilson GF linear-bend components fix the CH3CN fixed-conversion failure.
  Evidence: The Wilson GF validation path now appends two perpendicular linear-bend components for near-linear primitive bends without mutating the default Stage 3D primitive coordinate list. For CH3CN, selected indices change to `1;2;4;5;6;10;11;12;13;14;19;20`, where `19` and `20` are `linear_bend_component` rows for `C1-C2-N3`. `outputs\wilson_gf_batch_ch3cn\CH3CN__wilson_gf_validation.csv` reports `PASS`, max relative error `2.948152719245503e-07`, G-condition `33.4146494976297`, F-condition `53.48211609003476`, and warning `linear_bend_coordinate_used`.

- Observation: Aniline has fewer positive modes than the non-linear `3N-6` expected vibrational rank, even though the positive GF/ORCA mode counts match each other and the fixed-conversion comparison passes.
  Evidence: `outputs\wilson_gf_expand_aniline\aniline__wilson_gf_basis_diagnostics.csv` reports `expected_vibrational_rank` 36, `positive_orca_mode_count` 35, and `positive_gf_eigenvalue_count` 35. The validation CSV reports `PASS`, max relative error `5.6175518138189824e-08`, and warnings `positive_orca_mode_count_below_expected_vibrational_rank; positive_gf_eigenvalue_count_below_expected_vibrational_rank`.

- Observation: Acetanilide exposed a large-system conditioned-basis gap after the linear-bend and small-system conditioned-basis fixes.
  Evidence: The initial `outputs\wilson_gf_heavy_acetanilide\acetanilide__wilson_gf_validation.csv` run reported `FAIL`, basis size 51, expected rank 51, G-rank 50, G-condition `2349768983412.9233`, F-rank 51, F-condition `220982845170337.12`, positive ORCA modes 50, max relative error `6.054848166520849e-05`, and warnings `basis_rank_below_expected; g_ill_conditioned; f_ill_conditioned; positive_orca_mode_count_below_expected_vibrational_rank; positive_gf_eigenvalue_count_below_expected_vibrational_rank`. A singular-vector diagnostic showed the weakest mass-weighted B direction was dominated by three C5-centered bend rows: `ang(C4-C5-H9)`, `ang(C6-C5-H9)`, and `ang(C4-C5-C6)`.

- Observation: A deterministic mass-weighted residual-pivot validation basis fixes the acetanilide G-rank failure without loosening rank thresholds.
  Evidence: After adding the large-system fallback, `outputs\wilson_gf_heavy_acetanilide\acetanilide__wilson_gf_validation.csv` reports `PASS`, basis size 51, expected rank 51, G-rank 51, G-condition `1303.6329834917049`, F-rank 51, F-condition `5026.982264219395`, positive ORCA modes 50, max relative error `5.660912987811707e-08`, and warnings only for positive mode counts below expected rank.

- Observation: N-methylaniline exposed that G-only fallback acceptance can leave a severely ill-conditioned F matrix.
  Evidence: Before extending fallback acceptance, `outputs\wilson_gf_heavy_n_methylaniline\N-methylaniline__wilson_gf_validation.csv` reported `PASS`, G-rank 45, G-condition `578040400998.682`, F-rank 45, F-condition `64189060354395.125`, max relative error `6.788418311786068e-06`, and warning `f_ill_conditioned`. A mass-weighted pivot candidate gave G-condition `1326.4122537120506`, F-condition `33201.90313372605`, and max relative error `7.075733291133912e-08`.

- Observation: After F-condition-aware fallback acceptance, N-methylaniline no longer carries `f_ill_conditioned`.
  Evidence: The rerun `outputs\wilson_gf_heavy_n_methylaniline\N-methylaniline__wilson_gf_validation.csv` reports `PASS`, basis size 45, expected rank 45, G-rank 45, G-condition `1326.4122537120506`, F-rank 45, F-condition `33201.90313372605`, positive ORCA modes 44, max relative error `7.075733291133912e-08`, and warnings only for positive mode counts below expected rank.

- Observation: The remaining positive-mode-count warnings are caused by one negative ORCA vibrational mode in each affected `.hess`, not by missing parser rows.
  Evidence: `aniline.hess` has 42 frequencies, expected rank 36, six exact zero frequencies at indices 0-5, and one negative frequency at index 6 (`-367.3568575130425` cm-1). `N-methylaniline.hess` has 51 frequencies, expected rank 45, six exact zero frequencies at indices 0-5, and one negative frequency at index 6 (`-198.14705984767406` cm-1). `acetanilide.hess` has 57 frequencies, expected rank 51, six exact zero frequencies at indices 0-5, and one negative frequency at index 6 (`-49.58924532088554` cm-1).

- Observation: The closed GF validation mirrors the same one-mode nonpositive count in the affected cases.
  Evidence: The rerun basis diagnostics CSVs report `gf_nonpositive_eigenvalue_count` 1 for aniline, N-methylaniline, and acetanilide. The corresponding minimum nonpositive GF eigenvalues are `-0.01823750951305075`, `-0.00530596637314424`, and `-0.0003323259303203356`.

## Decision Log

- Decision: Implement the upgrade as a separate module, `src/wilson_gf.py`, instead of extending `src/ped.py` first.
  Rationale: The new work solves a GF eigenproblem and validates frequency reproduction; the current `src/ped.py` contains diagnostic PED projections. Keeping the prototype separate makes the numerical assumptions, unit conversion, and validation status easier to inspect and easier to disable.
  Date/Author: 2026-05-14 / Codex

- Decision: Gate pipeline integration behind a new opt-in CLI flag, `--wilson-gf-validation`.
  Rationale: Existing Stage 3D, PED v1, PED v2, Wilson-style PED audit, composed PED outputs, and viewer payloads must remain unchanged unless the user explicitly asks for the prototype outputs.
  Date/Author: 2026-05-14 / Codex

- Decision: Use H2O as the first numerical gate before validating NH3, formaldehyde, and ethene.
  Rationale: H2O has three atoms and three vibrational modes, making basis rank, mode mapping, and unit diagnostics easier to debug before larger molecules introduce mixed modes and more internal-coordinate choices.
  Date/Author: 2026-05-14 / Codex

- Decision: Do not promote wording to strict Wilson GF, full Wilson GF PED, VEDA-equivalent PED, or proven frequency reproduction until fixed unit conversion passes the target molecule set.
  Rationale: An empirical ratio between ORCA frequencies and `sqrt(lambda_GF)` is useful for diagnostics, but using ORCA frequencies to fit the conversion cannot by itself prove independent frequency reproduction.
  Date/Author: 2026-05-14 / Codex

- Decision: Preserve the ORCA normal-mode orientation rule everywhere the implementation touches normal modes.
  Rationale: ORCAVEDA stores normal modes as `normal_modes[:, mode]`; transposing this orientation would silently corrupt mode assignments and PED diagnostics.
  Date/Author: 2026-05-14 / Codex

- Decision: Failures must be reported as diagnostic rows or explicit warnings, not hidden behind empty `pass` blocks or broad exception suppression.
  Rationale: This is a scientific validation path. Singular bases, conversion failures, mapping mismatches, and ill-conditioned matrices are results that must remain visible.
  Date/Author: 2026-05-14 / Codex

- Decision: `src\wilson_gf.py` uses a mass-metric internal-force reconstruction for the closed GF eigenproblem instead of reusing `reconstruct_internal_force_matrix` from `src\ped.py`.
  Rationale: The existing PED helper is suitable for the Wilson GF-style diagnostic audit, but the closed GF eigenproblem requires the Wilson mass-metric back-transform `J = M^-1 B^T G^-1`. On H2O this reduces max relative error from about `2.30e-3` to about `1.82e-5`.
  Date/Author: 2026-05-14 / Codex

- Decision: Report `positive_orca_mode_count_below_expected_vibrational_rank` and `positive_gf_eigenvalue_count_below_expected_vibrational_rank` warnings when positive-mode counts are below the expected non-linear vibrational rank, without automatically changing `PASS` to `FAIL` if the compared positive ORCA and GF mode counts match and fixed conversion passes.
  Rationale: This preserves the meaning of the numerical validation actually performed while making incomplete-positive-mode cases explicit in machine-readable diagnostics.
  Date/Author: 2026-05-14 / Codex

- Decision: For large validation systems where exhaustive conditioned-basis search is infeasible, add a deterministic mass-weighted residual-pivot fallback that selects rows from the scaled mass-weighted Wilson B matrix and accept it only when it restores full G-rank and improves G-condition.
  Rationale: Acetanilide showed that the default selected basis can be full-rank in unscaled B but rank-deficient in G at the Wilson GF tolerance. The residual-pivot fallback directly targets the mass metric used by `G = B M^-1 B^T` and is restricted to the opt-in Wilson GF validation path.
  Date/Author: 2026-05-14 / Codex

- Decision: Extend large-system fallback acceptance to consider F-condition when the candidate basis keeps full G/F rank and improves the internal force matrix conditioning.
  Rationale: N-methylaniline had full G-rank and acceptable G-condition under the previous cutoff but retained an F-condition of `6.4189060354395125e13`. Since the closed GF validation depends on both `G` and the reconstructed internal `F`, the validation basis should not preserve a severely ill-conditioned F matrix when the deterministic mass-weighted pivot candidate improves it without rank loss.
  Date/Author: 2026-05-14 / Codex

- Decision: Keep positive-mode-count warnings as `PASS` diagnostics when positive ORCA and GF counts match and fixed conversion passes, but add explicit nonpositive-mode fields and warnings to the opt-in Wilson GF CSVs.
  Rationale: The affected source `.hess` files contain one negative vibrational frequency in addition to six zero translational/rotational modes. This is source evidence that the positive count is below `3N-6`; hiding it or converting it to a parser failure would be unsupported.
  Date/Author: 2026-05-14 / Codex

## Outcomes & Retrospective

Phase 0 and the standalone Phase 1 prototype are implemented. H2O basis formation is valid for the current non-linear target: basis indices `[0, 1, 2]`, B-rank `3`, B-condition `110.8498941006996`, G-rank `3`, G-condition `2.3108020332111234`, F-rank `3`, F-condition `16.098506422459817`.

H2O fixed-conversion validation is currently `PASS` for the standalone prototype. The fixed SI conversion produced max relative error `1.818626413830463e-05`, below the planned `1.0e-4` gate. The empirical ratio diagnostic is nearly constant (`median 2720.1791780974722`, `std 3.3457373459153256e-07`) and remains diagnostic-only rather than proof beyond the H2O gate.

Phase 3 opt-in integration is implemented. The CLI flag `--wilson-gf-validation` writes `H2O__wilson_gf_validation.csv`, `H2O__wilson_gf_ped_audit.csv`, and `H2O__wilson_gf_basis_diagnostics.csv` for the H2O CLI run under `outputs\wilson_gf_h2o_cli`. Default `analyze_orca_ped_like` runs remain without the new `wilson_gf_*` table keys.

Phase 4 target CLI runs completed on 2026-05-14. H2O, NH3, formaldehyde, and ethene pass fixed-conversion validation through the opt-in path. The ethene pass depends on the small-system Wilson-GF-conditioned basis fallback, so this remains validation-prototype evidence rather than a VEDA-equivalence claim.

| Molecule | Output directory | Basis size | Expected rank | G-rank | G-condition | F-rank | F-condition | Positive ORCA modes | Positive GF eigenvalues | Max relative error | Empirical ratio median | Empirical ratio std | Status | Warnings |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| H2O | `outputs\wilson_gf_h2o` | 3 | 3 | 3 | 2.3108020332111234 | 3 | 16.098506422459817 | 3 | 3 | 1.818626413830463e-05 | 2720.1791780974722 | 3.3457373459153256e-07 | PASS | none |
| NH3 | `outputs\wilson_gf_nh3` | 6 | 6 | 6 | 2.315274433059329 | 6 | 14.029291999970486 | 6 | 6 | 5.616351591931107e-08 | 2720.2284947182916 | 2.9370862370825758e-09 | PASS | none |
| formaldehyde | `outputs\wilson_gf_formaldehyde` | 6 | 6 | 6 | 1293895.778784718 | 6 | 754289.3497871537 | 6 | 6 | 5.666234451789424e-08 | 2720.2284947048374 | 5.033772184914088e-07 | PASS | none |
| ethene | `outputs\wilson_gf_ethene` | 12 | 12 | 12 | 28.556444875327514 | 12 | 65.64477191436387 | 12 | 12 | 5.61625148639572e-08 | 2720.228494729605 | 6.610274925267368e-09 | PASS | none |

Broader validation batch completed on 2026-05-14. This batch is evidence for the prototype's current numerical envelope, not for VEDA equivalence. After opt-in linear-bend augmentation for near-linear bends, all eight broader-batch molecules pass fixed-conversion validation; CH3CN carries the diagnostic warning `linear_bend_coordinate_used`.

| Molecule | Output directory | Basis size | Expected rank | G-rank | G-condition | F-rank | F-condition | Positive ORCA modes | Max relative error | Empirical ratio median | Empirical ratio std | Status | Warnings |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| MeOH | `outputs\wilson_gf_batch_meoh` | 12 | 12 | 12 | 27.457797487906166 | 12 | 236.09060650034456 | 12 | 5.616420455075335e-08 | 2720.2284947383223 | 5.294429644752303e-08 | PASS | none |
| EtOH | `outputs\wilson_gf_batch_etoh` | 21 | 21 | 21 | 68.05073106259707 | 21 | 251.2245599416959 | 21 | 5.616697269282427e-08 | 2720.2284947275293 | 2.587721576221904e-08 | PASS | none |
| Acetone | `outputs\wilson_gf_batch_acetone` | 24 | 24 | 24 | 586815398.5892428 | 24 | 8386033069.547947 | 24 | 5.625625863693999e-08 | 2720.228494754702 | 5.165432595206878e-05 | PASS | none |
| CH3CN | `outputs\wilson_gf_batch_ch3cn` | 12 | 12 | 12 | 33.4146494976297 | 12 | 53.48211609003476 | 12 | 2.948152719245503e-07 | 2720.2284947126723 | 0.0002652499396613242 | PASS | `linear_bend_coordinate_used` |
| DMSO | `outputs\wilson_gf_batch_dmso` | 24 | 24 | 24 | 73.13515964560295 | 24 | 109.88705956258443 | 24 | 5.651559409564404e-08 | 2720.228494727352 | 2.34840549448295e-07 | PASS | none |
| acetamide | `outputs\wilson_gf_batch_acetamide` | 21 | 21 | 21 | 286158581.4001269 | 21 | 5541794646.3023 | 20 | 5.6389704513833395e-08 | 2720.2284947332896 | 3.0801958065541683e-06 | PASS | none |
| phenol | `outputs\wilson_gf_batch_phenol` | 33 | 33 | 33 | 2108.3644479210197 | 33 | 5386.926054955553 | 33 | 5.6208317650343554e-08 | 2720.2284947158196 | 2.846752966781087e-08 | PASS | none |
| benzene | `outputs\wilson_gf_batch_benzene` | 30 | 30 | 30 | 2940.309649604426 | 30 | 3694.5604550950507 | 30 | 5.6171675848486615e-08 | 2720.2284947215594 | 1.027057300104927e-08 | PASS | none |

Larger-molecule validation expansion completed on 2026-05-14. This remains opt-in validation-prototype evidence only. All six runs report `PASS`, but several runs carry diagnostics that should be preserved in any future promotion decision.

| Molecule | Output directory | Basis size | Expected rank | G-rank | G-condition | F-rank | F-condition | Positive ORCA modes | Max relative error | Empirical ratio std | Status | Warnings |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| pyridine | `outputs\wilson_gf_expand_pyridine` | 27 | 27 | 27 | 890.9833031044437 | 27 | 1318.3128132775255 | 27 | 5.6182786600516276e-08 | 2.0781507999519033e-08 | PASS | none |
| aniline | `outputs\wilson_gf_expand_aniline` | 36 | 36 | 36 | 3566.3371357170154 | 36 | 10221.71749316018 | 35 | 5.6175518138189824e-08 | 3.2657277088738505e-08 | PASS | `positive_orca_mode_count_below_expected_vibrational_rank; positive_gf_eigenvalue_count_below_expected_vibrational_rank` |
| benzonitrile | `outputs\wilson_gf_expand_benzonitrile` | 33 | 33 | 33 | 23584295.879570205 | 33 | 150175039.7379429 | 33 | 1.8507632961046087e-07 | 0.00011260737307538036 | PASS | `near_linear_bend_coordinate` |
| nitrobenzene | `outputs\wilson_gf_expand_nitrobenzene` | 36 | 36 | 36 | 684.7607585627828 | 36 | 643.5661664372692 | 36 | 5.645657068730809e-08 | 1.524848128869999e-07 | PASS | none |
| benzoic acid | `outputs\wilson_gf_expand_benzoic_acid` | 39 | 39 | 39 | 6391195.43036815 | 39 | 5235683.357372124 | 39 | 6.102709090773594e-08 | 2.0889356718249092e-06 | PASS | none |
| acetophenone | `outputs\wilson_gf_expand_acetophenone` | 45 | 45 | 45 | 721.1511202210136 | 45 | 522.3656280535778 | 45 | 5.681210961188783e-08 | 2.68610955771519e-07 | PASS | none |

Heavy validation batch completed on 2026-05-14. This batch exposed the current prototype's first clear larger-molecule failure in acetanilide. After the large-system mass-weighted residual-pivot fallback, all eight heavy-batch molecules report `PASS`. Acetanilide still carries positive-mode-count diagnostics, so this is a validation-basis fix rather than a claim that the molecule has all `3N-6` positive reported modes.

| Molecule | Output directory | Basis size | Expected rank | G-rank | G-condition | F-rank | F-condition | Positive ORCA modes | Max relative error | Empirical ratio std | Status | Warnings |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| NMP | `outputs\wilson_gf_heavy_nmp` | 42 | 42 | 42 | 32693189.82267202 | 42 | 7338681.180009891 | 42 | 5.623115220542311e-08 | 1.4374315378575958e-07 | PASS | none |
| N-methylaniline | `outputs\wilson_gf_heavy_n_methylaniline` | 45 | 45 | 45 | 1326.4122537120506 | 45 | 33201.90313372605 | 44 | 7.075733291133912e-08 | 5.916851178679829e-06 | PASS | `positive_orca_mode_count_below_expected_vibrational_rank; positive_gf_eigenvalue_count_below_expected_vibrational_rank` |
| piperidine | `outputs\wilson_gf_heavy_piperidine` | 45 | 45 | 45 | 235.8380018026121 | 45 | 12.88624371545388 | 45 | 5.61720943171254e-08 | 1.221115751659388e-08 | PASS | none |
| cyclohexane chair | `outputs\wilson_gf_heavy_cyclohexane_chair` | 48 | 48 | 48 | 192.36761966245606 | 48 | 12.078000854648094 | 48 | 5.622735744573172e-08 | 2.7151408832149106e-08 | PASS | none |
| acetanilide | `outputs\wilson_gf_heavy_acetanilide` | 51 | 51 | 51 | 1303.6329834917049 | 51 | 5026.982264219395 | 50 | 5.660912987811707e-08 | 1.7287996782479167e-07 | PASS | `positive_orca_mode_count_below_expected_vibrational_rank; positive_gf_eigenvalue_count_below_expected_vibrational_rank` |
| monoethanolamine dimer cyclic | `outputs\wilson_gf_heavy_mea_dimer_cyclic` | 60 | 60 | 60 | 3630.5702406229834 | 60 | 4418.419103866282 | 60 | 5.815906969891593e-08 | 1.043531947424966e-06 | PASS | none |
| monoethanolamine dimer NH-to-O | `outputs\wilson_gf_heavy_mea_dimer_nh_to_o` | 60 | 60 | 60 | 64266.18861443776 | 60 | 159937.3061546995 | 60 | 6.347330456417468e-08 | 3.236394770993942e-06 | PASS | none |
| monoethanolamine dimer OH-to-N | `outputs\wilson_gf_heavy_mea_dimer_oh_to_n` | 60 | 60 | 60 | 349912.9918617817 | 60 | 998117.478527149 | 60 | 6.232690848679069e-08 | 5.435162527030724e-06 | PASS | none |

## Context and Orientation

`src/orca_parser.py` parses ORCA `.hess` files into `HessData`, including atoms, masses, Angstrom coordinates, frequencies, IR intensities, normal modes, and the Cartesian Hessian when the `$hessian` block is available. The ORCA block matrix parser must keep support for single-column headers.

`src/internal_coordinates.py` builds primitive and chemically annotated internal coordinates from atoms, coordinates, bonds, fragments, hydrogen bonds, and functional groups. These internal coordinates are the rows used to build the Wilson B matrix.

`src/b_matrix.py` computes the finite-difference B matrix and selects independent coordinate rows. The current basis selection is rank-oriented and PED-localization-aware; it is useful for diagnostics, but this plan must verify whether the selected basis is numerically suitable for GF diagonalization.

`src/ped.py` currently provides PED v1, PED v2, and Wilson GF-style PED diagnostics. PED v1 is a normalized B-matrix projection. PED v2 is a force-aware B/Hessian diagnostic. The Wilson-style audit builds `G`, reconstructs an internal force matrix from the Cartesian Hessian through a pseudoinverse, and reports projected potential-energy contributions for ORCA normal modes. These are existing evidence layers, not a closed GF frequency validation module.

`src/ORCAVEDA_patched_stage3D_v5_0.py` is the main pipeline entry point. It builds internal coordinates, B matrices, independent bases, Stage 3D assignment outputs, PED outputs, composed PED diagnostics, and report tables. `src/orcaveda_cli.py` owns argument parsing and must pass any new opt-in flag through to the pipeline without changing default behavior.

ORCAVEDA's current baseline remains Stage 3D v5.0 plus diagnostic PED layers. Stage 3D may be described as a geometric and weighted independent-coordinate assignment audit. It must not be renamed as strict VEDA PED or full Wilson GF PED.

## Plan of Work

Phase 0 is audit and prototype preparation. Confirm the existing basis selection, B matrix shape, Hessian availability, and current Wilson-style PED implementation. Record any mismatch between the external proposal and source evidence in this plan before coding. This phase ends when the implementer can state the exact basis size, rank, condition number, and available positive ORCA frequencies for H2O.

Phase 1 creates `src/wilson_gf.py` with core numerical functions and diagnostics. The module should expose a result dataclass, symmetric square-root decomposition for symmetric positive semidefinite matrices, a GF diagonalization function, a frequency-validation dataframe builder, a closed-PED diagnostic dataframe builder, and a basis-diagnostics dataframe builder if needed. The module must validate shapes, ranks, finite values, and condition numbers explicitly.

Phase 2 validates H2O before pipeline promotion. Build H2O internals using the same chemistry annotation path as the pipeline, compute B, select a `3N-6` basis, solve the GF eigenproblem, and compare sorted positive GF eigenvalues to positive ORCA vibrational frequencies. Record both a fixed physical conversion attempt and the empirical ratio diagnostic. If the fixed conversion fails but the empirical ratio is stable, the result is `WARN`, not `PASS`.

Phase 3 adds gated pipeline integration. Add `--wilson-gf-validation` to `src/orcaveda_cli.py`, pass it through to `run_orca_ped_like` and `analyze_general_hess_files`, and emit new CSV tables only when the flag is set. Existing output schemas and viewer payloads remain unchanged. If the prototype cannot run for a molecule, write diagnostic status rows where possible and continue producing existing ORCAVEDA outputs.

Phase 4 validates NH3, formaldehyde, and ethene. Run the same command path on all four target molecules and update this plan with observed max relative errors, conversion status, basis diagnostics, and PED normalization behavior. Only after these runs pass with fixed conversion should wording be promoted beyond "validation prototype."

## Concrete Steps

Work from the repository root:

    C:\Users\unive\Documents\Projects\orcaveda

Step 1: inspect the current implementation and record evidence in this plan.

    rg -n "build_wilson_g_matrix|reconstruct_internal_force_matrix|build_wilson_ped_audit_dataframe" src\ped.py
    rg -n "def select_independent_coordinates|def finite_difference_B" src\b_matrix.py
    rg -n "def cli_main|add_argument|parse_known_args" src\orcaveda_cli.py src\ORCAVEDA_patched_stage3D_v5_0.py

Expected: the current Wilson-style audit and CLI entry points are identified before code changes.

Step 2: create `src/wilson_gf.py`.

The module should include:

- `WILSON_GF_VALIDATION_METHOD`, with wording that says validation prototype and not VEDA-equivalent.
- `WilsonGFResult`, containing filename, atom count, Cartesian size, internal basis size, vibrational mode count, basis indices, G/F rank and condition, eigenvalues, eigenvectors, ORCA vibrational frequencies, reconstructed frequencies for fixed conversion, empirical ratio diagnostics, max relative errors, validation status, and warnings.
- `symmetric_sqrt_decomp(A, tol=1.0e-12)`.
- `wilson_gf_diagonalization(hess, internals, B, selected_idx, *, tol=1.0e-12, frequency_tol_relative=1.0e-4)`.
- `build_wilson_gf_validation_dataframe(result)`.
- `wilson_gf_closed_ped(result, hess, internals, B, selected_idx, *, top_n=8, tol=1.0e-12)`.
- `build_wilson_gf_basis_diagnostics_dataframe(result)` if the validation and PED tables do not already expose enough basis diagnostics.

Implementation requirements:

- Validate `hess.cartesian_hessian is not None`.
- Validate `B.shape[0] == len(internals)` and `B.shape[1] == 3 * len(hess.atoms)`.
- Use `expected_vibrational_rank = 3N-6` for current non-linear target molecules. Linear molecule handling may remain an explicit `TODO` warning because no target molecule in this plan is linear.
- Do not expect six zero eigenvalues from the internal-coordinate GF matrix.
- Reuse `wilson_coordinate_scales`, `build_wilson_g_matrix`, and `reconstruct_internal_force_matrix` from `src/ped.py` unless evidence shows their conventions are unsuitable.
- Solve the symmetric eigenproblem using `G^(1/2) F G^(1/2)`.
- Sort positive GF eigenvalues and positive ORCA vibrational frequencies consistently; record the mapping method in the output.
- Preserve `normal_modes[:, mode]` if normal-mode vectors are used anywhere in the implementation.
- Use explicit warnings such as `missing_cartesian_hessian`, `basis_rank_below_expected`, `g_ill_conditioned`, `fixed_conversion_failed`, `empirical_ratio_only`, and `mode_count_mismatch`.

Step 3: create `tests/test_wilson_gf.py`.

Focused tests should cover:

- `symmetric_sqrt_decomp` on identity and a known symmetric positive definite matrix.
- `symmetric_sqrt_decomp` rejection of non-square or non-symmetric input.
- Synthetic GF eigenvalues using known `G` and `F` matrices, without ORCAVEDA chemistry dependencies.
- H2O pipeline setup using the same annotation, internal-coordinate, B-matrix, and basis-selection path as `analyze_general_hess_files`.
- H2O validation dataframe contains method, basis diagnostics, conversion status, max relative error, and warnings columns.
- H2O closed-PED diagnostic rows have per-mode normalization sums close to 100% when valid eigenvectors are available. This proves implemented normalization, not VEDA equivalence.

Step 4: run syntax and focused tests.

    .\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py
    .\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q

Expected: syntax check passes. Focused tests either pass or fail with actionable numeric diagnostics recorded in this plan. Do not weaken tolerances without recording the observed error and rationale.

Step 5: integrate the opt-in CLI and pipeline path.

Add `--wilson-gf-validation` to `src/orcaveda_cli.py`, pass it to `run_orca_ped_like`, and thread it into `analyze_general_hess_files`. When disabled, the tables dictionary and manifest should match current behavior except for unrelated ordering only if unavoidable.

When enabled and `$hessian` is available, write:

- `{prefix}__wilson_gf_validation.csv`
- `{prefix}__wilson_gf_ped_audit.csv`
- `{prefix}__wilson_gf_basis_diagnostics.csv` if implemented as a separate table

If `$hessian` is missing or the basis is invalid, write a validation diagnostic row where feasible and continue existing outputs. Do not silently skip the prototype without a warning.

Step 6: run focused regression tests.

    .\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q

Expected: existing behavior remains unchanged when `--wilson-gf-validation` is absent.

Step 7: run target molecule validations.

    .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\H2O_freq.hess --outdir outputs\wilson_gf_h2o --wilson-gf-validation
    .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\NH3.hess --outdir outputs\wilson_gf_nh3 --wilson-gf-validation
    .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\formaldehyde.hess --outdir outputs\wilson_gf_formaldehyde --wilson-gf-validation
    .\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\ethene.hess --outdir outputs\wilson_gf_ethene --wilson-gf-validation

For each output directory, inspect the new CSV tables and record:

- basis size, rank, and condition number;
- number of positive ORCA vibrational modes;
- number of positive GF eigenvalues;
- fixed conversion max relative error;
- empirical ratio median and standard deviation;
- validation status (`PASS`, `WARN`, or `FAIL`);
- any warnings.

## Validation and Acceptance

The implementation is accepted as a Wilson GF validation prototype only when:

- `src/wilson_gf.py` compiles.
- `tests/test_wilson_gf.py` runs and records H2O diagnostics.
- H2O produces `wilson_gf_validation.csv` through the opt-in pipeline path.
- Existing focused PED and Stage 3D regression tests pass without the opt-in flag.
- With `--wilson-gf-validation`, all four target molecules produce validation CSV rows containing basis rank, basis condition, conversion method, mapping method, max relative error, empirical ratio diagnostics, and warnings.
- Existing CSV schemas are unchanged for `assignment_audit`, `ped_audit`, `ped_v2_force_audit`, `wilson_ped_audit`, composed PED outputs, `ped_stage3d_agreement`, and `ped_final_assignment`.
- No output, manifest, docstring, or report claims strict Wilson GF, full Wilson GF PED, VEDA equivalence, or proven frequency reproduction unless fixed physical conversion passes H2O, NH3, formaldehyde, and ethene with `max_relative_error < 1.0e-4`.

If fixed physical conversion does not pass but empirical ratios are stable, the prototype may be considered diagnostically useful with `WARN` status. That result must not be promoted to strict Wilson GF frequency reproduction.

## Idempotence and Recovery

Creating `src/wilson_gf.py`, `tests/test_wilson_gf.py`, and new opt-in tables is additive. Re-running the target molecule commands may overwrite files under `outputs/wilson_gf_*`, which is acceptable for generated validation artifacts.

The default ORCAVEDA run without `--wilson-gf-validation` is the recovery path. If the prototype fails for a molecule, existing Stage 3D, PED v1, PED v2, Wilson-style PED audit, composed PED, and viewer outputs must still be produced.

If a pipeline integration change causes existing tests to fail without the opt-in flag, revert only the integration path and keep the standalone prototype and focused tests for debugging.

If a numerical implementation fails H2O, do not broaden to NH3/formaldehyde/ethene until the failure is diagnosed. Record the exact observed eigenvalues, conversion diagnostics, basis rank/condition, and max relative error in `Outcomes & Retrospective`.

## Artifacts and Notes

Use this section to record implementation evidence as it is produced.

Initial source review placeholders:

- `src/ped.py`: Wilson-style diagnostic functions reviewed; exact line references to be added during implementation.
- `src/b_matrix.py`: independent-coordinate selection reviewed; exact H2O selected basis diagnostics to be added.
- `src/orcaveda_cli.py`: CLI threading reviewed; final flag behavior to be recorded.

Command transcript placeholders:

- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py`: passed on 2026-05-14.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q`: passed on 2026-05-14, latest run `8 passed in 24.57s`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q`: passed on 2026-05-14 after adding positive-mode-count warnings, latest run `9 passed in 24.68s`.
- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py tests\test_wilson_gf.py`: passed on 2026-05-14 after adding positive-mode-count warnings.
- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py tests\test_wilson_gf.py`: passed on 2026-05-14 after adding the large-system mass-weighted pivot fallback.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q`: passed on 2026-05-14 after adding the acetanilide large-system fallback test, latest run `10 passed in 28.76s`.
- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py tests\test_wilson_gf.py`: passed on 2026-05-14 after extending fallback acceptance to F-condition.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q`: passed on 2026-05-14 after adding the N-methylaniline F-condition fallback test, latest run `11 passed in 24.36s`.
- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py tests\test_wilson_gf.py`: passed on 2026-05-14 after adding nonpositive-mode diagnostics.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_wilson_gf.py -q`: passed on 2026-05-14 after adding nonpositive-mode diagnostic assertions, latest run `11 passed in 22.92s`.
- `.\.venv312\Scripts\python.exe -m py_compile src\wilson_gf.py src\orcaveda_cli.py src\ORCAVEDA_patched_stage3D_v5_0.py`: passed on 2026-05-14.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py -q`: passed on 2026-05-14, `19 passed in 1.31s`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q`: passed on 2026-05-14, latest run `21 passed in 28.08s`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q`: passed on 2026-05-14 after adding positive-mode-count warnings, latest run `21 passed in 28.37s`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q`: passed on 2026-05-14 after adding the large-system mass-weighted pivot fallback, latest run `21 passed in 26.31s`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q`: passed on 2026-05-14 after extending fallback acceptance to F-condition, latest run `21 passed in 25.93s`.
- `.\.venv312\Scripts\python.exe -m pytest tests\test_ped.py tests\test_stage3d_outputs.py tests\test_regression_baseline_outputs.py -q`: passed on 2026-05-14 after adding nonpositive-mode diagnostics, latest run `21 passed in 25.80s`.
- H2O CLI validation run: `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\H2O_freq.hess --outdir outputs\wilson_gf_h2o_cli --wilson-gf-validation`, completed on 2026-05-14.
- H2O Phase 4 validation run: `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\H2O_freq.hess --outdir outputs\wilson_gf_h2o --wilson-gf-validation`, completed on 2026-05-14.
- NH3 validation run: `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\NH3.hess --outdir outputs\wilson_gf_nh3 --wilson-gf-validation`, completed on 2026-05-14.
- Formaldehyde validation run: `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\formaldehyde.hess --outdir outputs\wilson_gf_formaldehyde --wilson-gf-validation`, completed on 2026-05-14.
- Ethene validation run: `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\ethene.hess --outdir outputs\wilson_gf_ethene --wilson-gf-validation`, completed on 2026-05-14.
- Broader validation batch runs completed on 2026-05-14 for `MeOH_freq.hess`, `EtOH_freq.hess`, `Acetone_freq.hess`, `CH3CN_freq.hess`, `DMSO_freq.hess`, `acetamide.hess`, `phenol.hess`, and `benzene.hess`, each through `src\ORCAVEDA_patched_stage3D_v5_0.py --wilson-gf-validation`.
- Larger-molecule validation expansion runs completed on 2026-05-14 for `pyridine.hess`, `aniline.hess`, `benzonitrile.hess`, `nitrobenzene.hess`, `benzoic_acid.hess`, and `acetophenone.hess`, each through `src\ORCAVEDA_patched_stage3D_v5_0.py --wilson-gf-validation`.
- Heavy validation batch runs completed on 2026-05-14 for `NMP_freq.hess`, `N-methylaniline.hess`, `piperidine.hess`, `cyclohexane_chair.hess`, `acetanilide.hess`, `monoethanolamine_dimer_cyclic_DFT.hess`, `monoethanolamine_dimer_NH_to_O_DFT.hess`, and `monoethanolamine_dimer_OH_to_N_DFT.hess`, each through `src\ORCAVEDA_patched_stage3D_v5_0.py --wilson-gf-validation`.
- Larger-molecule and heavy validation batches were rerun on 2026-05-14 after adding the mass-weighted pivot fallback. All rerun CLI commands exited 0.
- Larger-molecule and heavy validation batches were rerun on 2026-05-14 after extending fallback acceptance to F-condition. All rerun CLI commands exited 0.
- Affected warning-case validation outputs were rerun on 2026-05-14 after adding nonpositive-mode diagnostics for `aniline.hess`, `N-methylaniline.hess`, and `acetanilide.hess`; all rerun CLI commands exited 0.

CSV evidence placeholders:

- H2O `wilson_gf_validation.csv`: generated by the opt-in CLI run as `outputs\wilson_gf_h2o_cli\H2O__wilson_gf_validation.csv`; all validation rows report `PASS`, max relative error `0.000018`, no warnings.
- H2O `wilson_gf_ped_audit.csv`: generated by the opt-in CLI run as `outputs\wilson_gf_h2o_cli\H2O__wilson_gf_ped_audit.csv`.
- H2O `wilson_gf_basis_diagnostics.csv`: generated by the opt-in CLI run as `outputs\wilson_gf_h2o_cli\H2O__wilson_gf_basis_diagnostics.csv`.
- H2O Phase 4 CSVs: generated under `outputs\wilson_gf_h2o`; validation status `PASS`, max relative error `1.818626413830463e-05`, warnings none.
- NH3 CSVs: generated under `outputs\wilson_gf_nh3`; validation status `PASS`, max relative error `5.616351591931107e-08`, warnings none.
- Formaldehyde CSVs: generated under `outputs\wilson_gf_formaldehyde`; validation status `PASS`, max relative error `5.666234451789424e-08`, warnings none.
- Ethene CSVs: generated under `outputs\wilson_gf_ethene`; validation status `PASS`, selected indices `0;4;5;6;7;9;12;13;14;15;17;18`, max relative error `5.61625148639572e-08`, warnings none.
- Broader batch CSVs: generated under `outputs\wilson_gf_batch_*`; MeOH, EtOH, acetone, DMSO, acetamide, phenol, benzene, and CH3CN report `PASS`; CH3CN carries `linear_bend_coordinate_used`.
- Larger-molecule expansion CSVs after the mass-weighted pivot fallback: generated under `outputs\wilson_gf_expand_*`; pyridine, aniline, benzonitrile, nitrobenzene, benzoic acid, and acetophenone report `PASS`. Warnings remain for aniline positive-mode counts below expected rank and benzonitrile near-linear bend.
- Heavy batch CSVs after F-condition-aware fallback acceptance: generated under `outputs\wilson_gf_heavy_*`; NMP, N-methylaniline, piperidine, cyclohexane chair, acetanilide, and the three monoethanolamine dimers report `PASS`. N-methylaniline and acetanilide retain positive-mode-count diagnostics below expected rank; no heavy-batch molecule currently reports `basis_rank_below_expected`, `g_ill_conditioned`, or `f_ill_conditioned`.
- Nonpositive-mode diagnostics: `outputs\wilson_gf_expand_aniline\aniline__wilson_gf_basis_diagnostics.csv`, `outputs\wilson_gf_heavy_n_methylaniline\N-methylaniline__wilson_gf_basis_diagnostics.csv`, and `outputs\wilson_gf_heavy_acetanilide\acetanilide__wilson_gf_basis_diagnostics.csv` now include `orca_nonpositive_mode_count`, `orca_min_nonpositive_frequency_cm-1`, `gf_nonpositive_eigenvalue_count`, and `gf_min_nonpositive_eigenvalue`. The affected rows carry `nonpositive_orca_modes_within_expected_vibrational_space` and `nonpositive_gf_eigenvalues_within_expected_vibrational_space`.
- CH3CN exhaustive primitive-basis diagnostic: evaluated 1136 full-rank 12-coordinate subsets; passing subsets `0`; best max relative error `0.01845477680255841`.
- CH3CN linear-bend CLI rerun: `.\.venv312\Scripts\python.exe src\ORCAVEDA_patched_stage3D_v5_0.py data\hess\CH3CN_freq.hess --outdir outputs\wilson_gf_batch_ch3cn --wilson-gf-validation`, completed on 2026-05-14 with `PASS`, max relative error `2.948152719245503e-07`.

Standalone H2O diagnostics from `src\wilson_gf.py`:

- basis indices: `[0, 1, 2]`
- B-rank / condition: `3` / `110.8498941006996`
- G-rank / condition: `3` / `2.3108020332111234`
- F-rank / condition: `3` / `16.098506422459817`
- GF eigenvalues: `0.320797023521, 2.21605760856, 2.31699381851`
- ORCA positive frequencies: `1540.680823, 4049.375352, 4140.568268`
- reconstructed fixed-conversion frequencies: `1540.708842, 4049.448994, 4140.643569`
- max relative error: `1.818626413830463e-05`
- empirical ratio median / std: `2720.1791780974722` / `3.3457373459153256e-07`
- status / warnings: `PASS` / none

## Interfaces and Dependencies

New public module:

- `src/wilson_gf.py`

Proposed public functions and data structures:

- `WilsonGFResult`
- `symmetric_sqrt_decomp(A, tol=1.0e-12) -> tuple[np.ndarray, np.ndarray]`
- `wilson_gf_diagonalization(hess, internals, B, selected_idx, *, tol=1.0e-12, frequency_tol_relative=1.0e-4) -> WilsonGFResult`
- `build_wilson_gf_validation_dataframe(result) -> pandas.DataFrame`
- `wilson_gf_closed_ped(result, hess, internals, B, selected_idx, *, top_n=8, tol=1.0e-12) -> pandas.DataFrame`
- `build_wilson_gf_basis_diagnostics_dataframe(result) -> pandas.DataFrame`, optional if diagnostics are not already sufficient in the validation table

New CLI flag:

- `--wilson-gf-validation`

New optional output schemas:

`{prefix}__wilson_gf_validation.csv` must include at minimum:

- `Source`
- `Filename`
- `mode_index`
- `orca_frequency_cm-1`
- `gf_eigenvalue`
- `reconstructed_frequency_cm-1`
- `fixed_conversion_relative_error`
- `empirical_ratio_frequency_cm1_per_sqrt_lambda`
- `mapping_method`
- `conversion_method`
- `validation_status`
- `max_relative_error`
- `empirical_ratio_median`
- `empirical_ratio_std`
- `orca_nonpositive_mode_count`
- `orca_min_nonpositive_frequency_cm-1`
- `gf_nonpositive_eigenvalue_count`
- `gf_min_nonpositive_eigenvalue`
- `basis_size`
- `expected_vibrational_rank`
- `g_rank`
- `g_condition`
- `f_rank`
- `f_condition`
- `warnings`
- `method`

`{prefix}__wilson_gf_ped_audit.csv` must include at minimum:

- `Source`
- `Filename`
- `mode`
- `frequency_cm-1`
- `gf_rank`
- `coord_index`
- `internal_coordinate`
- `coordinate_kind`
- `coordinate_family`
- `signed_ped_fraction`
- `contribution_percent`
- `normalization_sum_percent`
- `basis_size`
- `validation_status`
- `max_relative_error`
- `warnings`
- `method`

`{prefix}__wilson_gf_basis_diagnostics.csv`, if emitted separately, must include:

- `Source`
- `Filename`
- `basis_size`
- `expected_vibrational_rank`
- `selected_indices`
- `g_rank`
- `g_condition`
- `f_rank`
- `f_condition`
- `positive_orca_mode_count`
- `positive_gf_eigenvalue_count`
- `orca_nonpositive_mode_count`
- `orca_min_nonpositive_frequency_cm-1`
- `gf_nonpositive_eigenvalue_count`
- `gf_min_nonpositive_eigenvalue`
- `warnings`

No new dependency is allowed beyond NumPy, Pandas, and existing ORCAVEDA code unless this plan is updated with a specific rationale and validation impact.

Existing interfaces that must remain unchanged by default:

- `assignment_audit.csv`
- `ped_audit.csv`
- `ped_v2_force_audit.csv`
- `wilson_ped_audit.csv`
- `composed_ped_audit.csv`
- `composed_ped_v2_force_audit.csv`
- `composed_wilson_ped_audit.csv`
- `ped_stage3d_agreement.csv`
- `ped_final_assignment.csv`
- interactive viewer HTML/JSON payload fields
- `expectations/regression_expectations_stage3D_v5_0.json`
