# Changelog

All notable changes to this reproducibility release are documented here.
This project adheres to [Semantic Versioning](https://semver.org/) for its
analysis pipeline and to the manuscript's own revision numbering.

## [1.0.1] - 2026-10-28

Revision additions requested during the JNS submission polish (supplementary
naming and inline cross-references; the manuscript now cites the supplement as
`Supplementary Material` rather than `Online Resource 1`).

### Added
- `code/44_age_sex_and_threshold_sensitivity.R` — R `logistf` Firth fits for
  Supplementary Table S8 (age/sex-adjusted resection models: Taiwan Knosp OR
  4.08 [2.19–9.08], Spain 7.53 [3.70–19.37]) and for the Supplementary Table S4
  lower block (binary LOS threshold >7 d Taiwan / >4 d Spain, adjusted for
  previous surgery and for age+sex). Engine, seed (20260925), penalized
  profile-likelihood CI and penalized LR p-method identical to
  `13_canonical_firth.R`.
- `results/sr_age_sex_threshold_sensitivity.csv` — reference values for the
  published Knosp term of the six revision models above (see the file's
  `status` and `note` columns). Committed as a check-value snapshot; re-running
  `44_*` against the raw public files regenerates the full coefficient table.
- `CHANGELOG.md` (this file).

### Changed
- `code/00_prep_public.py` now propagates a harmonised `male` column (1 = male,
  0 = female) into `data/d3_clean.csv` so the Taiwan age/sex Firth model is
  reproducible from the released pipeline. Sex / gender is detected across
  common raw header variants (Sex, Gender, M/F, Gen); if none is present the
  script warns and leaves `male` as `NaN`, in which case script `44` skips the
  Taiwan age/sex fit.
- `README.md` — added a "Revision analyses (post-v1.0.0)" section, listed script
  `44` in the repository layout and inserted `Rscript
  code/44_age_sex_and_threshold_sensitivity.R` in the "How to run" pipeline.
- `CITATION.cff` — version bumped to 1.0.1, date-released 2026-10-28.

### Not changed
- No statistical model, cleaning rule, estimand, seed, or figure logic for the
  v1.0.0 analyses was altered. The France (D1) cohort remains removed.

## [1.0.0] - 2026-10-28

Initial reproducibility release for the two-cohort JNS secondary analysis.

### Added
- Full analysis pipeline: `00_prep_public.py`, `07_d3_functionalform.py`,
  `08_d4_prep.py`, `09_d4_domain_separation.py`, `13_canonical_firth.R`,
  `d3_firth_logistf.R`, `d4_firth_logistf.R`, `firth_utils.py`,
  `33_pooled_nb2_and_mi_sensitivity.py`, `38_endpoint_concordance_and_loo.py`,
  `40_fig1_two_cohort.py`, `41_adjusted_absolute_risk_v2.py`,
  `42_figures_v4.py`, `43_verify_from_raw.py`.
- Canonical aggregate outputs under `results/`, Figures 1–3 (PNG + PDF) and the
  LOS leave-one-out supplementary panel under `figures/`.
- `data/DATA_MANIFEST.md` documenting the two raw public source files and their
  previous-surgery codings.
- `environment/` (pinned requirements + session_info).
