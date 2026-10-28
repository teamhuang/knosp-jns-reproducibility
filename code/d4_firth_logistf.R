# d4_firth_logistf.R — External validation of the anatomical->resectability domain in
# D4 (Mendeley NFPA multicentre cohort), using the SAME Firth specification as D3:
# residual (GTR Parcial) ~ Knosp + previous surgery, with penalised profile-likelihood
# CIs and penalised likelihood-ratio P values (logistf).
suppressMessages(library(logistf))
ROOT <- Sys.getenv("KNOSP_REPO_ROOT", ".")
RES  <- file.path(ROOT, "results")
d4 <- read.csv(file.path(ROOT, "data", "d4_clean.csv"))

fit_tbl <- function(fit, lab, n){
  co <- fit$coefficients
  df <- data.frame(model=lab, n=n, term=names(co), OR=exp(as.numeric(co)),
                   OR_lo=exp(as.numeric(fit$ci.lower)), OR_hi=exp(as.numeric(fit$ci.upper)),
                   p=as.numeric(fit$prob), ci_method="profile penalized-likelihood",
                   p_method="penalized likelihood-ratio", stringsAsFactors=FALSE)
  df[df$term != "(Intercept)", ]
}
run <- function(form, lab){
  f <- suppressMessages(logistf(form, data=d4, profilelik=TRUE))
  fit_tbl(f, lab, nrow(model.frame(form, d4)))
}

out <- list()
out[[1]] <- run(residual ~ knosp + prev_surgery, "D4_VALID_Firth_ModelA")
out[[2]] <- run(residual ~ knosp_high + prev_surgery, "D4_Firth_threshold34")
out[[3]] <- run(residual ~ knosp + prev_surgery + log2vol, "D4_Firth_ModelB_voladj")
d4m <- do.call(rbind, out)
write.csv(d4m, file.path(RES, "d4_firth_logistf_models.csv"), row.names=FALSE)
cat("=== D4 external validation (Firth, residual ~ Knosp + prev_surgery) ===\n")
print(d4m[d4m$model=="D4_VALID_Firth_ModelA", c("term","OR","OR_lo","OR_hi","p")], row.names=FALSE)
cat("\n=== all D4 models ===\n")
print(d4m[, c("model","n","term","OR","OR_lo","OR_hi","p")], row.names=FALSE)
