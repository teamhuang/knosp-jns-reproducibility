# d3_firth_logistf.R — Authoritative Firth inference for D3 via logistf.
# Firth OR + penalized PROFILE-likelihood CI (f$ci.lower/upper) + penalized
# likelihood-ratio P (f$prob). FLIC calibrated probabilities via f$flic (Puhr 2017).
suppressMessages(library(logistf))
# ROOT defaults to the repo root; run from the repo directory, or set
#   KNOSP_REPO_ROOT=/path/to/repo  and run  Rscript code/d3_firth_logistf.R
ROOT <- Sys.getenv("KNOSP_REPO_ROOT", ".")
RES  <- file.path(ROOT, "results")

d3 <- read.csv(file.path(ROOT, "data", "d3_clean.csv"))
d3$Size <- tolower(trimws(as.character(d3$Size)))
d3$is_macro <- d3$Size %in% c("macro", "giant", "huge")
d3$Knosp_high <- as.integer(d3$Knosp >= 3)
d3$residual <- as.integer(d3$residual); d3$subtotal_eor <- as.integer(d3$subtotal_eor)

fit_tbl <- function(fit, lab, outcome, cohort, n, ev){
  co <- fit$coefficients
  df <- data.frame(model=lab, outcome=outcome, cohort=cohort, n=n, events=ev,
                   term=names(co), beta=as.numeric(co), OR=exp(as.numeric(co)),
                   OR_lo=exp(as.numeric(fit$ci.lower)), OR_hi=exp(as.numeric(fit$ci.upper)),
                   p=as.numeric(fit$prob), ci_method="profile penalized-likelihood",
                   p_method="penalized likelihood-ratio", stringsAsFactors=FALSE)
  df[df$term != "(Intercept)", ]
}
run_model <- function(dat, form, lab, outcome, cohort){
  f <- tryCatch(suppressMessages(logistf(form, data=dat, profilelik=TRUE)), error=function(e) NULL)
  if (is.null(f)) return(NULL)
  fit_tbl(f, lab, outcome, cohort, nrow(dat), sum(model.frame(form, dat)[,1]))
}
analytic <- function(dat){ dat[complete.cases(dat[, c("Knosp","prev_surgery","log2vol","residual","subtotal_eor")]), ] }

out <- list()
for (y in c("residual","subtotal_eor")){
  full <- analytic(d3); macro <- full[full$is_macro, ]
  out[[length(out)+1]] <- run_model(full,  as.formula(paste(y,"~ Knosp + prev_surgery")), "PRIMARY_Firth_ModelA_full", y, "full_n100")
  out[[length(out)+1]] <- run_model(macro, as.formula(paste(y,"~ Knosp + prev_surgery")), "Firth_ModelA_macro", y, "macro_n80")
  out[[length(out)+1]] <- run_model(full,  as.formula(paste(y,"~ Knosp + prev_surgery + log2vol")), "Firth_ModelB_full", y, "full_n100")
  out[[length(out)+1]] <- run_model(macro, as.formula(paste(y,"~ Knosp + prev_surgery + log2vol")), "Firth_ModelB_macro", y, "macro_n80")
  out[[length(out)+1]] <- run_model(full,  as.formula(paste(y,"~ Knosp_high + prev_surgery")), "Firth_threshold_full", y, "full_n100")
  gm <- glm(as.formula(paste(y,"~ Knosp + prev_surgery")), data=full, family=binomial)
  cs <- summary(gm)$coefficients
  out[[length(out)+1]] <- data.frame(model="MLE_ModelA_full", outcome=y, cohort="full_n100", n=nrow(full),
                    events=sum(full[[y]]), term=rownames(cs)[-1], beta=cs[-1,1], OR=exp(cs[-1,1]),
                    OR_lo=exp(cs[-1,1]-1.96*cs[-1,2]), OR_hi=exp(cs[-1,1]+1.96*cs[-1,2]),
                    p=cs[-1,4], ci_method="Wald", p_method="Wald", stringsAsFactors=FALSE)
}
d3_models <- do.call(rbind, out)
write.csv(d3_models, file.path(RES, "d3_firth_logistf_models.csv"), row.names=FALSE)
cat("=== D3 Firth (logistf) models: residual ===\n")
print(d3_models[d3_models$outcome=="residual", c("model","cohort","term","OR","OR_lo","OR_hi","p")], row.names=FALSE)

# ---- functional form: nested Firth penalized-LR (linear vs categorical Knosp) ----
ff <- data.frame()
for (y in c("residual","subtotal_eor")){
  full <- analytic(d3)
  for (coh in c("full_n100","macro_n80")){
    dat <- if (coh=="full_n100") full else full[full$is_macro,]
    f_lin <- suppressMessages(logistf(as.formula(paste(y,"~ Knosp + prev_surgery")), data=dat))
    f_cat <- suppressMessages(logistf(as.formula(paste(y,"~ factor(Knosp) + prev_surgery")), data=dat))
    lr <- 2*(f_cat$loglik[1] - f_lin$loglik[1]); df <- 3
    p <- max(0, pchisq(lr, df, lower.tail=FALSE))
    ff <- rbind(ff, data.frame(outcome=y, cohort=coh, n=nrow(dat), events=sum(dat[[y]]),
                               LR_linear_vs_cat=round(lr,3), df=df, p_nonlinearity=round(p,4),
                               test="penalized likelihood-ratio (nested Firth)", stringsAsFactors=FALSE))
  }
}
write.csv(ff, file.path(RES, "d3_knosp_functional_form.csv"), row.names=FALSE)
cat("\n=== Knosp functional-form (nested Firth penalized LR) ===\n"); print(ff, row.names=FALSE)

# ---- FLIC calibrated probabilities per Knosp grade (primary full Model A) ----
# FLIC (Puhr 2017): keep Firth slopes, correct the intercept so the mean predicted
# probability equals the observed event rate. lp from the fitted Firth model.
full <- analytic(d3)
f <- suppressMessages(logistf(residual ~ Knosp + prev_surgery, data=full))
lp <- as.numeric(f$linear.predictors)
rate <- mean(full$residual)
a <- uniroot(function(k) mean(1/(1+exp(-(lp+k)))) - rate, c(-10,10))$root
flic_p <- 1/(1+exp(-(lp+a)))
flp <- data.frame(knosp=full$Knosp, flic_prob=flic_p, firth_prob=1/(1+exp(-lp)), residual=full$residual)
agg <- aggregate(cbind(flic_prob, firth_prob) ~ knosp, data=flp, FUN=mean)
obsr <- aggregate(residual ~ knosp, data=flp, FUN=mean)
cnt  <- aggregate(residual ~ knosp, data=flp, FUN=length)
agg$observed_rate <- obsr$residual; agg$n <- cnt$residual
write.csv(agg, file.path(RES, "d3_flic_probs.csv"), row.names=FALSE)
cat("\n=== FLIC probabilities (primary full Model A) ===\n"); print(agg, row.names=FALSE)
cat("\nintercept correction a=", round(a,3), " observed overall rate=", round(rate,3),
    " FLIC mean=", round(mean(flic_p),3), "\n")
cat("\n[logistf D3 done]\n")
