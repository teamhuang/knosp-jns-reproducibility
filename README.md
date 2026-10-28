# Knosp grade and incomplete resection — cross-cohort reproducibility (JNS)

Reproducible analysis code and canonical model outputs for:

> **Cross-Cohort Reproducibility of Knosp Grade Associations With Incomplete
> Resection After Pituitary Surgery: A Secondary Analysis of Two Independent
> Cohorts.** *Target journal: Journal of Neurosurgery (JNS).*

This is a **secondary analysis of two independently collected, publicly
available, de-identified individual-level pituitary-surgery cohorts** (Taiwan =
index; Spain = replication). It tests the **reproducibility of associations**,
not the development or external validation of a prediction model, and is
reported per **STROBE**. No new patient data were generated or recruited.

---

## Headline claims this repository backs

| Claim | Reported value | Reproduced from raw public files? |
|---|---|---|
| Knosp → incomplete resection, Taiwan (Firth OR / 95% CI) | 3.97 (2.22–8.40) | point estimate ✓ (<0.1%); CI = R `logistf` canonical |
| Knosp → incomplete resection, Spain (Firth OR / 95% CI) | 5.69 (3.10–12.21) | point estimate ✓; CI = R `logistf` canonical |
| Pooled common-slope Knosp OR | 5.04 (3.23–8.55) | point estimate ✓; CI = R `logistf` canonical |
| Pooled Knosp × cohort interaction OR | 1.49 (0.57–3.93) | point estimate ✓; CI = R `logistf` canonical |
| Model-adjusted absolute risk difference, Taiwan (Knosp 3–4 vs 0–2) | 42.8 pp (26.8–59.9) | ✓ 1000-bootstrap |
| Model-adjusted absolute risk difference, Spain (Knosp 3–4 vs 0–2) | 59.7 pp (46.0–73.0) | ✓ 1000-bootstrap |
| Knosp → continuous LOS (NB2 IRR), Taiwan / Spain | 1.24 (1.12–1.38) / 0.98 (0.90–1.08) | ✓ |

Every non-CI quantity is re-derived from the raw public files by
`code/43_verify_from_raw.py` (see *Verification* below). Firth confidence
intervals are taken from the **canonical R `logistf` engine** (penalized
profile likelihood); the Python Firth routine reproduces the point estimates to
<0.1% but its profile-likelihood intervals are systematically narrower and are
**not** used for any reported interval — this is recorded per-row in
`results/verification_manifest_v4.csv`.

---

## Repository layout

```
.
├── code/                     # analysis pipeline (Python 3.10 + R 4.6)
│   ├── _paths.py             # repo-relative path resolution (env-overridable)
│   ├── 00_prep_public.py     # Taiwan (d3) analytic table from the raw PLOS workbook
│   ├── 07_d3_functionalform.py       # Knosp functional form, threshold, collinearity, bootstrap (Supp.)
│   ├── 08_d4_prep.py         # Spain (d4) analytic table from the raw Mendeley .sav
│   ├── 09_d4_domain_separation.py    # within-cohort anatomy→resection vs anatomy→LOS
│   ├── 13_canonical_firth.R  # CANONICAL Firth models → firth_{d3,d4,pooled}_logistf.csv
│   ├── d3_firth_logistf.R    # Taiwan sensitivity/FLIC detail (logistf)
│   ├── d4_firth_logistf.R    # Spain Firth detail (logistf)
│   ├── firth_utils.py        # Python Firth MLE + penalized profile CI (cross-check engine)
│   ├── 33_pooled_nb2_and_mi_sensitivity.py   # pooled Knosp×cohort NB2 interaction sensitivity
│   ├── 38_endpoint_concordance_and_loo.py    # endpoint concordance + leave-one-out influence
│   ├── 40_fig1_two_cohort.py # Figure 1 schematic
│   ├── 41_adjusted_absolute_risk_v2.py       # corrected g-computation adjusted absolute risk
│   ├── 42_figures_v4.py      # Figures 1–3 from canonical result tables
│   ├── 43_verify_from_raw.py # CLEAN-ROOM verifier: reproduces headline numbers from raw
│   └── 44_age_sex_and_threshold_sensitivity.R  # Supp. Tables S4 threshold + S8 age/sex Firth (revision)
├── results/                  # canonical AGGREGATE model outputs (safe to publish; no patient rows)
├── figures/                  # Figure1–3 (PNG 300dpi + vector PDF) and figS6
├── data/                     # DATA_MANIFEST.md only (raw files are public + separately licensed)
├── environment/              # requirements.txt + session_info.md (pinned versions)
└── LICENSE / CITATION.cff / README.md
```

## Why the raw data are not included

Both cohorts are **third-party patient-level clinical data** released by their
original authors under their own licences. Redistributing individual-level rows
here is neither necessary nor appropriate. Instead, the exact source files,
DOIs, licences and column dictionary are documented in
[`data/DATA_MANIFEST.md`](data/DATA_MANIFEST.md). Place the two files below in
`data/` (or point the environment variables at them) and the pipeline runs
end-to-end:

- `data/PONE-D-21-14669R1-Raw_data.xlsx` — Taiwan (Chen P.-Y., figshare,
  doi:10.6084/m9.figshare.15042771)
- `data/D4_mendeley_nfpa.sav` — Spain (Pérez-Lopez C., Mendeley Data,
  doi:10.17632/kb55bvtv3y.1, CC BY 4.0)

## Environment

- **Python 3.10**, dependencies pinned in [`environment/requirements.txt`](environment/requirements.txt)
  (`pandas`, `numpy`, `scipy`, `statsmodels`, `pyreadstat`, `openpyxl`, `matplotlib`).
- **R 4.6** with **`logistf` 1.26.1** for the canonical Firth fits (penalized
  profile-likelihood CIs and penalized LR p-values). See
  [`environment/session_info.md`](environment/session_info.md).
- **Fixed seed `20260925`** for every bootstrap/randomisation.

## How to run

From the repository root (paths default to `./data`, `./results`, `./figures`;
override with `KNOSP_DATA_DIR`, `KNOSP_RESULTS_DIR`, `KNOSP_FIG_DIR`):

```bash
# 0) put the two public raw files in ./data  (see data/DATA_MANIFEST.md)

# 1) build the two analytic tables (patient-level; written locally, git-ignored)
python code/00_prep_public.py
python code/08_d4_prep.py

# 2) CANONICAL Firth models (authoritative intervals) — requires R + logistf
Rscript code/13_canonical_firth.R

# 3) exploratory / sensitivity tables
python code/07_d3_functionalform.py
python code/09_d4_domain_separation.py
python code/33_pooled_nb2_and_mi_sensitivity.py
python code/38_endpoint_concordance_and_loo.py
python code/41_adjusted_absolute_risk_v2.py
Rscript code/44_age_sex_and_threshold_sensitivity.R

# 4) figures
python code/40_fig1_two_cohort.py
python code/42_figures_v4.py

# 5) independent clean-room verification from the raw public files only
python code/43_verify_from_raw.py
```

## Verification (what we ran before release)

Running `code/43_verify_from_raw.py` **directly from a fresh checkout** (reading
only the two public raw files) reproduces **26/26 non-confidence-interval
headline quantities**, with maximum absolute difference **0.05**, and regenerates
`results/verification_manifest_v4.csv` and
`results/macroadenoma_selection_audit.csv` **byte-identically** to the committed
files. Scripts `07`, `33`, `38`, `41` likewise regenerate their committed
outputs identically. The verifier prints:

```
Non-CI quantities reproduced from raw: 26/26 pass | max abs_diff = 0.0519
CI rows flagged False by design (canonical engine = R logistf, not Python approximation).
```

The 10 "False" rows are exactly the Firth **intervals** and are expected:
`results/verification_manifest_v4.csv` labels their canonical engine as
`R logistf (canonical)`; the Python Firth profile CI is an approximation and is
never used for a reported interval.

## Scope note — removed cohort and portability changes

- **French cohort removed at revision.** An earlier three-cohort draft included a
  French database (internal id `d1`) and a French multiple-imputation LOS
  sensitivity (`sr_mi_linearity_sensitivity.json`). That cohort was **deleted
  entirely**; no French code path, data, or result is part of this release.
  Where a working script also built the French branch, the release copy was
  **re-scoped to the two retained cohorts** with the statistical logic of those
  two cohorts left unchanged (this is documented inline in `00_prep_public.py`
  and `33_pooled_nb2_and_mi_sensitivity.py`).
- **Path portability only.** The working scripts originally hard-coded a single
  machine's absolute paths; the release resolves paths repo-relative via
  `_paths.py` (env-overridable). No model, cleaning rule, estimand, seed, or
  figure logic was altered — `41` additionally reads the derived Spain
  previous-surgery codings (`prev_grouped/prev_strict/prev_ambl2no`) emitted by
  `08`, matching the values behind the shipped models.
- **Adjusted absolute risk uses the corrected estimand** (individual-level
  standardization `E[expit(Xβ)]`, previous surgery standardized to cohort
  prevalence), not the expit of a band mean — `code/41_adjusted_absolute_risk_v2.py`.

## Revision analyses (post-v1.0.0)

Two supplementary-table analyses were added after the initial v1.0.0 tag and are
included in v1.0.1:

- **`code/44_age_sex_and_threshold_sensitivity.R`** — (A) Firth models of
  incomplete resection additionally adjusted for age (per 10-year increase) and
  sex (manuscript Supplementary Table S8; Taiwan OR 4.08, Spain OR 7.53), and
  (B) Firth models of binary LOS thresholds (>7 d Taiwan, >4 d Spain) adjusted
  for previous surgery and for age+sex (manuscript Supplementary Table S4 lower
  block). Engine, seed, and CI/P method are identical to the canonical
  `13_canonical_firth.R` fits.
- **`code/00_prep_public.py`** now propagates a harmonised `male` (1=male,
  0=female) column into `data/d3_clean.csv`, detecting the sex/gender header
  across common raw variants (Sex, Gender, M/F, Gen). If none is found the
  script warns and leaves `male` as `NaN`; script `44` then skips the Taiwan
  age/sex fit and prints the reason.

Because the two public cohorts are third-party and are **not redistributed
here**, `results/sr_age_sex_threshold_sensitivity.csv` committed with v1.0.1
carries only the **published Knosp-term values** quoted in the manuscript as a
reference output; running `code/44_age_sex_and_threshold_sensitivity.R` against
the raw files (after `00_prep_public.py` and `08_d4_prep.py`) regenerates the
full coefficient table and overwrites the reference file. Users should treat
the reference rows as a check-value snapshot, not as a substitute for the
authoritative `logistf` re-fit.

## Figures

`figures/` contains the three main figures (300 dpi PNG + editable vector PDF).
Panel/forest values are read from the canonical tables in `results/`
(no hard-coded effect sizes). `figures/supplementary/figS6_los_leave_one_out.png`
is the leave-one-out influence diagnostic.

## Licence & citation

Code and aggregate results are released under **CC BY 4.0** to match the most
restrictive upstream data licence; the original cohort data remain under their
own licences (see `data/DATA_MANIFEST.md`). Please cite the manuscript and the
source cohort datasets. See [`CITATION.cff`](CITATION.cff).

## Persistent DOI

This repository is intended to be archived as a **versioned release with a
version-specific DOI** (e.g. via the GitHub–Zenodo integration, which mints a new
DOI per tagged release rather than a concept DOI). Once the release is tagged and
archived, paste the resolved version DOI into the manuscript's **Code
availability** statement.
