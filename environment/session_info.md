# Software environment

## Python (used for preparation, NB2 LOS, bootstrap, verification, figures)

- Python **3.10**
- numpy 1.26.4 · pandas 2.2.2 · scipy 1.14.1 · statsmodels 0.14.6 ·
  pyreadstat 1.3.6 · openpyxl 3.1 · matplotlib 3.8.4

Fixed random seed for all bootstrap / randomisation steps: **20260925**.
Numpy's PCG64 generator (`np.random.default_rng(20260925)`) is used for the
patient bootstrap and the g-computation standardisation.

## R (CANONICAL Firth engine)

- R **4.6**
- **`logistf` 1.26.1** — Firth penalized logistic regression with
  **penalized profile-likelihood confidence intervals** and **penalized
  likelihood-ratio p-values** (Puhr 2017). These are the intervals and p-values
  reported in the manuscript (Tables and forest plots); `code/13_canonical_firth.R`,
  `code/d3_firth_logistf.R`, `code/d4_firth_logistf.R`.

## Engine note (reimplementation boundary)

The Python Firth implementation in `code/firth_utils.py` is provided as an
independent **cross-check**. It reproduces the R `logistf` **point estimates to
<0.1%**, but its penalized profile-likelihood intervals are **systematically
narrower** and are therefore **not** used for any reported interval. This is
flagged per-row in `results/verification_manifest_v4.csv`
(`engine = "R logistf (canonical)"` vs `"Python re-run from raw"`).

Negative-binomial (NB2, joint MLE) LOS models use `statsmodels.NegativeBinomial`
(`loglike_method="nb2"`); robustness fits use HC3 OLS on log(LOS) and Poisson with
HC1 covariance, matching the manuscript Methods.
