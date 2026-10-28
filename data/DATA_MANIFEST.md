# Data manifest (raw public inputs — not redistributed here)

Both cohorts are **third-party, publicly available, de-identified individual-level
clinical datasets** released by their original authors. This repository contains
only the **code** and the **aggregate model outputs**; the raw files must be
downloaded from their sources and placed in `data/` (or pointed to via the
`KNOSP_TAIWAN_XLSX` / `KNOSP_SPAIN_SAV` environment variables) before running the
pipeline.

## 1. Taiwan — index cohort

| | |
|---|---|
| File | `data/PONE-D-21-14669R1-Raw_data.xlsx` |
| Source | Chen P.-Y. *PONE-D-21-14669R1-Raw_data.xlsx.* figshare; 2021. |
| DOI | **doi:10.6084/m9.figshare.15042771** |
| Structure | Two worksheets: `One hand` and `Two Hand` (surgical approach). |
| Original ethics | Chang Gung Memorial Hospital IRB **201901259B0**. |
| Raw n → analytic n | 103 → 100 (3 rows dropped for missing Knosp/previous-surgery/volume/residual). |

Key columns used (original header → harmonised name):

```
Previous OP   -> prev_surgery   (binary; previous pituitary surgery)
Tumor_volume  -> volume_cm3     (cm^3; log2vol = log2(volume))
MRI Residual  -> residual       (binary outcome: residual on post-op MRI = incomplete)
Hospital stay -> los_days       (length of stay, days; exploratory outcome)
Knosp         -> Knosp          (0-4 ordinal invasion grade; primary exposure)
EOR           -> subtotal_eor / incomplete_eor   (extent-of-resection category, alternate endpoints)
Diameter, Size                  (tumour maximum dimension (mm) and size category: micro/macro/giant/huge)
approach      -> one_hand / two_hand (derived from worksheet)
```

## 2. Spain — replication cohort

| | |
|---|---|
| File | `data/D4_mendeley_nfpa.sav` |
| Source | Pérez-Lopez C. *Data for: Isointense T2-weighted signal negatively influences the chances of gross total resection of nonfunctioning pituitary adenomas…* Mendeley Data; 2020. |
| DOI | **doi:10.17632/kb55bvtv3y.1** |
| Licence | **CC BY 4.0** (use complies with this licence). |
| Original ethics | Local ethics-committee approval with patient consent. |
| n | 100 non-functioning pituitary adenomas (single operating centre; multi-centre referral). |

Key columns used (raw variable → harmonised name), built by `code/08_d4_prep.py`:

```
RxKnosPre        -> knosp           (0-4)
GTR              -> gtr / residual  (GTR Completa=1 -> residual = 1 - gtr = incomplete)
QxEstanciaDias   -> los_days        (post-operative hospital length of stay; ICU stay is separate)
IQPrevia_Grupo   -> prev_grouped (= prev_surgery)  authors' pre-grouped 0/1 previous surgery (20 yes)
IQ_hipofisis_previa -> prev_strict / prev_ambl2no  coding sensitivity (codes {1=pituitary, 2/3=ambiguous})
RxVolPre_ml      -> volume_ml       (log2vol)
Hospital         -> referring institution (NOT the operating centre)
```

### Previous-surgery coding variants (Spain)

`IQ_hipofisis_previa` uses codes `{0=none, 1=pituitary surgery, 2/3=ambiguous other}`.
Three codings are used across the models:

| Variant | Definition | n | yes |
|---|---|---|---|
| `prev_grouped` (= authors' `IQPrevia_Grupo`) | codes 1/2/3 → yes (PRIMARY) | 100 | 20 |
| `prev_strict` | documented codes only; ambiguous 2/3 dropped | 95 | 15 |
| `prev_ambl2no` | ambiguous codes 2/3 recoded to No | 100 | 15 |

## Outcome-definition caveat (recorded in the manuscript)

The two cohorts operationalise "incomplete resection" differently — the Taiwan
cohort via post-operative **MRI residual**, the Spain cohort via **volumetric
GTR**. The Knosp direction replicates across these differing endpoints, but the
absolute levels and the two endpoints are **not** interchangeable; this is why the
paper frames the finding as a *reproducible association*, not a *transportable
absolute risk*. A per-grade endpoint-concordance table for the Taiwan cohort is
produced by `code/38_endpoint_concordance_and_loo.py`
(`results/sr_taiwan_endpoint_concordance.csv`).

## Macroadenoma sub-cohort reconciliation (Supplementary S6)

`results/macroadenoma_selection_audit.csv` records, row-by-row, the macroadenoma
counts by surgical approach derived two independent ways from the public Taiwan
workbook (Size category and maximum-dimension > 10 mm), which agree exactly
(80 macro/giant/huge of the 100 analysed), and flags where they differ from the
source publication's reported per-approach figure.

## Provenance / removed cohort

An earlier draft analysed a **third (French) cohort** (internal id `d1`,
`Database.csv`). That cohort was removed in full at revision. No French dataset,
variable, code path, or result appears in this release.
