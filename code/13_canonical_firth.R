# 13_canonical_firth.R — SINGLE canonical Firth engine (logistf) for the revised SR manuscript.
# Every Firth point estimate, penalised profile-likelihood 95% CI and penalised
# likelihood-ratio P value in the paper is produced here, so no model appears with
# two different intervals. Python is used only for non-Firth models (NB2/OLS/Poisson)
# and as an independent point-estimate cross-check.
suppressMessages(library(logistf))
suppressWarnings(Sys.setlocale("LC_ALL", "English"))
ROOT <- "."
RES  <- file.path(ROOT, "results")
d3 <- read.csv(file.path(ROOT,"data","d3_clean.csv"))
d4 <- read.csv(file.path(ROOT,"data","d4_clean.csv"))

tbl <- function(fit, lab, n){
  co <- fit$coefficients
  data.frame(model=lab, n=n, term=names(co),
             OR=exp(as.numeric(co)),
             OR_lo=exp(as.numeric(fit$ci.lower)),
             OR_hi=exp(as.numeric(fit$ci.upper)),
             p=as.numeric(fit$prob),
             ci_method="penalized profile-likelihood (logistf)",
             p_method="penalized likelihood-ratio (logistf)",
             stringsAsFactors=FALSE)
}
run <- function(form, dat, lab){ f <- logistf(form, data=dat, profilelik=TRUE); tbl(f, lab, nrow(model.frame(form, dat))) }

# ---- D3 (Taiwan index) ----
a3 <- subset(d3, !is.na(Knosp) & !is.na(prev_surgery) & !is.na(residual) & !is.na(log2vol))
out3 <- list(
  run(residual ~ Knosp + prev_surgery, a3, "D3_ModelA_residual"),
  run(residual ~ Knosp_high + prev_surgery, a3, "D3_threshold_residual"),
  run(residual ~ Knosp + prev_surgery + log2vol, a3, "D3_ModelB_residual"),
  run(subtotal_eor ~ Knosp + prev_surgery, a3, "D3_ModelA_subtotal"),
  run(residual ~ Knosp + prev_surgery, subset(a3, macro80==1), "D3_ModelA_macro80")
)
d3all <- do.call(rbind, out3); write.csv(d3all, file.path(RES,"firth_d3_logistf.csv"), row.names=FALSE)

# ---- D4 (Spanish replication) primary + coding-rule sensitivities ----
out4 <- list(
  run(residual ~ knosp + prev_grouped, d4, "D4_ModelA_grouped_n100"),
  run(residual ~ knosp, d4, "D4_KnospOnly_n100"),
  run(residual ~ knosp_high + prev_grouped, d4, "D4_threshold_n100"),
  run(residual ~ knosp + prev_grouped + log2vol, d4, "D4_ModelB_n100"),
  run(residual ~ knosp + prev_strict, subset(d4, !is.na(prev_strict)), "D4_sens_documented_n95"),
  run(residual ~ knosp + prev_ambl2no, d4, "D4_sens_ambiguousToNo_n100")
)
d4all <- do.call(rbind, out4); write.csv(d4all, file.path(RES,"firth_d4_logistf.csv"), row.names=FALSE)

# ---- Pooled (Taiwan + Spain) common-slope vs interaction ----
# Harmonised with the per-cohort primary Model A: adjust for Knosp + previous surgery + cohort.
tw <- data.frame(knosp=d3$Knosp, prev=d3$prev_surgery, resid=d3$residual, cohort="Taiwan")
es <- data.frame(knosp=d4$knosp,  prev=d4$prev_grouped, resid=d4$residual, cohort="Spain")
pool <- rbind(tw, es)
pool <- pool[!is.na(pool$knosp) & !is.na(pool$resid) & !is.na(pool$prev), ]
pool$cohort <- factor(pool$cohort, levels=c("Taiwan","Spain"))
pool$knosp <- as.numeric(pool$knosp); pool$prev <- as.numeric(pool$prev); pool$resid <- as.numeric(pool$resid)
common <- run(resid ~ knosp + prev + cohort, pool, "POOLED_common_slope")          # knosp = pooled common-slope OR (adj prev surgery + cohort)
inter  <- run(resid ~ knosp + prev + cohort + knosp:cohort, pool, "POOLED_interaction")  # knosp = Taiwan slope; knosp:cohortSpain = difference in slope
pooled <- rbind(common, inter)
write.csv(pooled, file.path(RES,"firth_pooled_logistf.csv"), row.names=FALSE)

cat("=== D3 (Taiwan) Firth ===\n"); print(d3all[d3all$term!="Intercept", c("model","n","term","OR","OR_lo","OR_hi","p")], row.names=FALSE)
cat("\n=== D4 (Spain) Firth + coding rules ===\n"); print(d4all[d4all$term!="Intercept", c("model","n","term","OR","OR_lo","OR_hi","p")], row.names=FALSE)
cat("\n=== Pooled common-slope vs interaction ===\n"); print(pooled[pooled$term!="Intercept", c("model","n","term","OR","OR_lo","OR_hi","p")], row.names=FALSE)
cat("\n[13 done] wrote firth_d3_logistf.csv, firth_d4_logistf.csv, firth_pooled_logistf.csv\n")
