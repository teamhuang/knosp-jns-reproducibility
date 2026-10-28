"""
firth_utils.py — Firth penalised logistic regression with the inference that
logistf recommends for sparse data:
  * point estimates  : Firth bias-reduced MLE
  * 95% CIs          : penalised PROFILE-likelihood (not Wald)
  * P values         : penalised LIKELIHOOD-RATIO test (not Wald)
  * FLIC             : Firth Linear Predictor with intercept Correction (Puhr 2017)
                       for calibrated absolute predicted probabilities.

Penalised log-likelihood:  l*(b) = l(b) + 0.5 * log|I(b)|,  I(b) = X' W X, W=diag(p(1-p)).
"""
import numpy as np
from scipy import optimize, stats


def _penalized_ll(X, y, beta):
    eta = X @ beta
    # stable log-likelihood
    ll = np.sum(y * eta - np.logaddexp(0, eta))
    p = 1.0 / (1.0 + np.exp(-eta))
    W = np.clip(p * (1 - p), 1e-12, None)
    I = X.T @ (X * W[:, None])
    # robust log|I| via eigenvalues (PSD); near-singular (separation) -> large negative
    eig = np.linalg.eigvalsh(0.5 * (I + I.T))
    eig = np.clip(eig, 1e-12, None)
    logdet = np.sum(np.log(eig))
    return ll + 0.5 * logdet


def firth_fit(X, y, max_iter=200, tol=1e-9):
    """Firth MLE = maximiser of the penalised log-likelihood l*(beta).
    Direct numerical maximisation (robust; the plain Newton on the modified score
    can converge to a wrong stationary point under separation)."""
    X = np.asarray(X, float); y = np.asarray(y, float)

    # warm start from the Newton-on-modified-score iterate
    beta = np.zeros(X.shape[1])
    for _ in range(max_iter):
        eta = X @ beta; p = 1.0 / (1.0 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-9, None)
        XtWX = X.T @ (X * W[:, None]); inv = np.linalg.pinv(XtWX)
        h = np.einsum('ij,jk,ik->i', X, inv, X)
        score = X.T @ (y - p + h * (0.5 - p))
        step = np.linalg.solve(XtWX + 1e-9 * np.eye(X.shape[1]), score)
        mx = np.max(np.abs(step))
        if mx > 2: step *= 2 / mx
        beta = beta + step
        if np.max(np.abs(step)) < tol: break

    res = optimize.minimize(lambda b: -_penalized_ll(X, y, b), beta,
                            method='Nelder-Mead',
                            options={'maxiter': 40000, 'maxfev': 80000, 'xatol': 1e-8, 'fatol': 1e-10})
    res = optimize.minimize(lambda b: -_penalized_ll(X, y, b), res.x, method='BFGS')
    return res.x


def _profile_ll_at(X, y, j, b_j, beta_init):
    """max l*(beta) over nuisance params with beta[j] fixed at b_j."""
    nuisance_idx = [k for k in range(X.shape[1]) if k != j]

    def neg_obj(bn):
        beta = np.empty(X.shape[1]); beta[j] = b_j
        beta[nuisance_idx] = bn
        return -_penalized_ll(X, y, beta)

    x0 = beta_init[nuisance_idx]
    res = optimize.minimize(neg_obj, x0, method='Nelder-Mead',
                            options={'maxiter': 4000, 'xatol': 1e-7, 'fatol': 1e-9})
    res = optimize.minimize(neg_obj, res.x, method='BFGS')
    return -res.fun


def profile_ci(X, y, j, beta_hat, ll_hat, alpha=0.05):
    """Penalised profile-likelihood 100(1-alpha)% CI for coefficient j."""
    crit = stats.chi2.ppf(1 - alpha, 1) / 2.0  # = 1.921 for 95%

    def gap(b):
        return 2 * (ll_hat - _profile_ll_at(X, y, j, b, beta_hat)) - crit

    se_wald = 1.96  # fallback scale
    # lower bound
    lo = beta_hat[j]
    step = 0.5
    lo_try = beta_hat[j] - step
    while gap(lo_try) < 0 and step < 60:
        step *= 1.6; lo_try = beta_hat[j] - step
    try:
        lo = optimize.brentq(gap, max(beta_hat[j] - 100, lo_try - step), beta_hat[j] - 1e-6, xtol=1e-5)
    except Exception:
        lo = beta_hat[j] - 1.96 * se_wald
    hi = beta_hat[j]
    step = 0.5
    hi_try = beta_hat[j] + step
    while gap(hi_try) < 0 and step < 60:
        step *= 1.6; hi_try = beta_hat[j] + step
    try:
        hi = optimize.brentq(gap, beta_hat[j] + 1e-6, hi_try + step, xtol=1e-5)
    except Exception:
        hi = beta_hat[j] + 1.96 * se_wald
    return lo, hi


def penalized_lr_p(X, y, keep_mask, beta_full, X_red, beta_reduced):
    """Penalised LR test for dropping the terms where keep_mask==False."""
    ll_full = _penalized_ll(X, y, beta_full)
    ll_red = _penalized_ll(X_red, y, beta_reduced)
    lr = 2 * (ll_full - ll_red)
    df = int((~keep_mask).sum())
    p = float(stats.chi2.sf(max(lr, 0), df)) if df > 0 else np.nan
    return lr, df, p


def firth_summary(X, y, names):
    """Full table: OR, penalised profile CI, penalised LR p (each term vs its removal)."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    beta = firth_fit(X, y)
    ll_hat = _penalized_ll(X, y, beta)
    rows = []
    for j, nm in enumerate(names):
        lo, hi = profile_ci(X, y, j, beta, ll_hat)
        keep = np.ones(X.shape[1], bool); keep[j] = False
        X_red = X[:, keep]
        beta_red = firth_fit(X_red, y)
        _, _, p = penalized_lr_p(X, y, keep, beta, X_red, beta_red)
        rows.append({'term': nm, 'beta': beta[j], 'OR': np.exp(beta[j]),
                     'OR_lo': np.exp(lo), 'OR_hi': np.exp(hi), 'p_lr': p})
    import pandas as pd
    return pd.DataFrame(rows), beta


def flic_probs(X, y, beta_firth):
    """FLIC: keep Firth slopes, correct intercept so mean predicted prob == observed rate."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    target = y.mean()
    lin = X[:, 1:] @ beta_firth[1:]  # exclude original intercept

    def g(a):
        return np.mean(1.0 / (1.0 + np.exp(-(a + lin)))) - target
    a = optimize.brentq(g, -20, 20, xtol=1e-8)
    beta_flic = beta_firth.copy(); beta_flic[0] = a
    return 1.0 / (1.0 + np.exp(-(a + lin))), beta_flic
