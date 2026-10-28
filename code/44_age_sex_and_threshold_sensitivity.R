# code/44_age_sex_and_threshold_sensitivity.R
# Revision analyses added post-v1.0.0: (A) age/sex-adjusted Firth models of
# incomplete resection (manuscript Supplementary Table S8) and (B) Firth models
# of binary LOS thresholds adjusted for previous surgery and for age+sex
# (manuscript Supplementary Table S4, threshold rows).
#
# Engine: R `logistf` 1.26.1, penalized profile-likelihood 95% CI, penalized
# likelihood-ratio P. Fixed seed 20260925 (logistf does not use RNG here but the
# seed is recorded for pipeline consistency).
#
# Requires code/00_prep_public.py (Taiwan d3_clean.csv, now propagates `male`)
# and code/08_d4_prep.py (Spain d4_clean.csv) to have been run first.
#
# Run from the repo root, or set KNOSP_REPO_ROOT=/path/to/repo.
suppressMessages(library(logistf))

ROOT <- Sys.getenv("KNOSP_REPO_ROOT", ".")
RES  <- file.path(ROOT, "results")
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE)
set.seed(20260925)

# ---- shared fit-table builder (mirrors d3_firth_logistf.R) -------------------
fit_tbl <- function(fit, lab, outcome, cohort, n, ev){
  co  <- fit$coefficients
  df  <- data.frame(model = lab, outcome = outcome, cohort = cohort,
                    n = n, events = ev,
                    term  = names(co),
                    beta  = as.numeric(co),
                    OR    = exp(as.numeric(co)),
                    OR_lo = exp(as.numeric(fit$ci.lower)),
                    OR_hi = exp(as.numeric(fit$ci.upper)),
                    p     = as.numeric(fit$prob),
                    ci_method = "profile penalized-likelihood",
                    p_method  = "penalized likelihood-ratio",
                    engine    = "R logistf 1.26.1",
                    seed      = 20260925L,
                    stringsAsFactors = FALSE)
  df[df$term != "(Intercept)", ]
}

run_or_null <- function(dat, form, lab, outcome, cohort){
  f <- tryCatch(suppressMessages(logistf(form, data = dat, profilelik = TRUE)),
                error = function(e) {
                  cat(sprintf("  skip %s [%s]: %s\n", lab, cohort, conditionMessage(e)))
                  NULL })
  if (is.null(f)) return(NULL)
  ev <- sum(model.frame(form, dat)[, 1], na.rm = TRUE)
  fit_tbl(f, lab, outcome, cohort, nrow(dat), ev)
}

out <- list()

# =============================================================================
# (A) Age/sex-adjusted Firth models of incomplete resection (Supp. Table S8)
# =============================================================================
d3 <- read.csv(file.path(ROOT, "data", "d3_clean.csv"), stringsAsFactors = FALSE)
d4 <- read.csv(file.path(ROOT, "data", "d4_clean.csv"), stringsAsFactors = FALSE)
d3$age10  <- d3$Age / 10
d4$age10  <- d4$age / 10
d4$Knosp  <- d4$knosp
d4$prev_surgery <- d4$prev_surgery
d4$residual <- d4$residual

# Taiwan: complete cases on Knosp, prev_surgery, Age, male, residual
d3_ok <- d3[stats::complete.cases(d3[, c("Knosp","prev_surgery","age10","male","residual")]), ]
if (nrow(d3_ok) > 0 && !all(is.na(d3_ok$male))){
  out[[length(out)+1]] <- run_or_null(
    d3_ok, residual ~ Knosp + prev_surgery + age10 + male,
    "AGESEX_Firth_ModelA_full", "MRI residual (Taiwan)", "full_n100")
} else {
  cat("  WARN: Taiwan `male` unavailable in d3_clean.csv; skipping D3 age/sex model.\n")
  cat("        Re-run code/00_prep_public.py after patching to propagate sex.\n")
}
# Spain
d4_ok <- d4[stats::complete.cases(d4[, c("knosp","prev_surgery","age10","male","residual")]), ]
out[[length(out)+1]] <- run_or_null(
  d4_ok, residual ~ knosp + prev_surgery + age10 + male,
  "AGESEX_Firth_ModelA_full", "incomplete volumetric resection (Spain)", "full_n100")

# =============================================================================
# (B) Threshold LOS Firth models (Supp. Table S4 lower block)
# =============================================================================
# Taiwan LOS > 7 days (post hoc)
d3_los <- d3[stats::complete.cases(d3[, c("Knosp","prev_surgery","los_gt7")]), ]
d3_los$los_gt7 <- as.integer(d3_los$los_gt7)
out[[length(out)+1]] <- run_or_null(
  d3_los, los_gt7 ~ Knosp + prev_surgery,
  "LOS_THRESHOLD_prevsurgery", "Taiwan LOS>7d", "full_n100")
d3_los_as <- d3[stats::complete.cases(d3[, c("Knosp","prev_surgery","los_gt7","age10","male")]), ]
d3_los_as$los_gt7 <- as.integer(d3_los_as$los_gt7)
if (nrow(d3_los_as) > 0 && !all(is.na(d3_los_as$male)))
  out[[length(out)+1]] <- run_or_null(
    d3_los_as, los_gt7 ~ Knosp + prev_surgery + age10 + male,
    "LOS_THRESHOLD_age_sex", "Taiwan LOS>7d", "full_n100")

# Spain LOS > 4 days (post hoc)
d4_los <- d4[stats::complete.cases(d4[, c("knosp","prev_surgery","plos_gt4")]), ]
d4_los$plos_gt4 <- as.integer(d4_los$plos_gt4)
out[[length(out)+1]] <- run_or_null(
  d4_los, plos_gt4 ~ knosp + prev_surgery,
  "LOS_THRESHOLD_prevsurgery", "Spain LOS>4d", "full_n100")
d4_los_as <- d4[stats::complete.cases(d4[, c("knosp","prev_surgery","plos_gt4","age10","male")]), ]
d4_los_as$plos_gt4 <- as.integer(d4_los_as$plos_gt4)
out[[length(out)+1]] <- run_or_null(
  d4_los_as, plos_gt4 ~ knosp + prev_surgery + age10 + male,
  "LOS_THRESHOLD_age_sex", "Spain LOS>4d", "full_n100")

# =============================================================================
res <- do.call(rbind, Filter(Negate(is.null), out))
rownames(res) <- NULL
write.csv(res, file.path(RES, "sr_age_sex_threshold_sensitivity.csv"), row.names = FALSE)
cat("\n=== wrote results/sr_age_sex_threshold_sensitivity.csv (", nrow(res), " rows ) ===\n")
print(res[, c("model","cohort","outcome","term","OR","OR_lo","OR_hi","p")], row.names = FALSE)
