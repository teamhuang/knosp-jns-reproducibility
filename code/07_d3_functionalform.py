"""
07_d3_descriptives.py — D3 descriptives + collinearity + Firth bootstrap (rev4).

Firth MODELS, functional-form LR and FLIC probabilities are produced by the
authoritative R/logistf script (d3_firth_logistf.R). This script only produces:
  * observed residual/subtotal rate per Knosp grade (full n=100 + macro n=80)
  * macroadenoma cohort reconciliation (our 80 vs source 78)
  * Knosp-volume correlation + VIF (estimand framing, not deletion)
  * bootstrap coefficient stability using the CORRECTED Firth estimator
"""
from pathlib import Path
import sys
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).parent))
from firth_utils import firth_fit

from _paths import DATA_DIR as DATA, RESULTS_DIR as RES
DATA = Path(DATA); RES = Path(RES)
d3 = pd.read_csv(DATA/'d3_clean.csv')
d3['size_lab'] = d3['Size'].astype(str).str.strip().str.lower()
d3['is_macro'] = d3['size_lab'].isin(['macro', 'giant', 'huge'])

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k/n; den = 1+z*z/n; c = (p+z*z/(2*n))/den
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return (max(0, c-h), min(1, c+h))

def frame(y, cohort):
    cols = ['Knosp','prev_surgery','log2vol',y]
    df = d3.dropna(subset=cols).copy()
    if cohort == 'macro':
        df = df[df['is_macro']]
    df[y] = df[y].astype(int)
    return df

# ---------- MODEL-SPECIFIC missingness audit (rev5) ----------
# Analytic n per model is set by the variables THAT model needs, not one blanket
# complete-case. Audit the raw 103 directly.
audit_rows = []
for lab, cols in [('PRIMARY Model A (MRI residual)', ['residual','Knosp','prev_surgery']),
                  ('SECONDARY Model A (Subtotal EOR)', ['subtotal_eor','Knosp','prev_surgery']),
                  ('Model B (volume-adjusted)', ['residual','Knosp','prev_surgery','log2vol'])]:
    m = d3[cols].notna().all(axis=1)
    audit_rows.append({'model':lab, 'raw_n':len(d3), 'analytic_n':int(m.sum()), 'dropped':int((~m).sum()),
                       'dropped_reason':'3 entirely-blank one-hand rows (ID/Size/Knosp/prev_surgery/volume/EOR all missing)'})
audit = pd.DataFrame(audit_rows)
audit.to_csv(RES/'d3_model_specific_missingness.csv', index=False)
print('=== model-specific missingness audit ===\n', audit.to_string(index=False))

# ---------- cohort reconciliation ----------
recon = pd.DataFrame([{'definition': 'Size in {macro,giant,huge} (max diameter >10mm)',
                       'our_n': int(d3['is_macro'].sum()),
                       'our_one_hand': int((d3['is_macro'] & (d3['approach']=='one_hand')).sum()),
                       'our_two_hand': int((d3['is_macro'] & (d3['approach']=='two_hand')).sum()),
                       'source_paper_n': 78, 'source_one_hand': 30, 'source_two_hand': 48,
                       'note': 'two-hand matches exactly; one-hand +2 not reproducible from public vars; '
                               'therefore full n=100 is PRIMARY, macro n=80 is sensitivity'}])
recon.to_csv(RES/'d3_macroadenoma_reconciliation.csv', index=False)
print('=== macroadenoma reconciliation ===\n', recon.to_string(index=False))

# ---------- observed rate per Knosp grade (both cohorts) ----------
rows = []
for y in ['residual','subtotal_eor']:
    for coh in ['full_n100','macro_n80']:
        df = frame(y, 'macro' if coh=='macro' else 'full')
        for g in [0,1,2,3,4]:
            sub = df[df['Knosp']==g]; k=int(sub[y].sum()); n=len(sub); lo,hi = wilson(k,n)
            rows.append({'outcome':y,'cohort':coh,'knosp_grade':g,'n':n,'events':k,
                         'rate':round(k/n,3) if n else np.nan,'ci_lo':round(lo,3),'ci_hi':round(hi,3)})
obs = pd.DataFrame(rows)
obs.to_csv(RES/'d3_knosp_observed_rates.csv', index=False)
print('\n=== observed rate by Knosp grade (full_n100 residual) ===')
print(obs[(obs.outcome=='residual')&(obs.cohort=='full_n100')].to_string(index=False))

# ---------- collinearity / VIF (full primary cohort) ----------
df = frame('residual','full')
sp = stats.spearmanr(df['Knosp'], df['log2vol'])
def vif(df, cols):
    X = sm.add_constant(df[cols].astype(float).to_numpy(), has_constant='add'); out={}
    for i,c in enumerate(cols):
        yj=X[:,i+1]; Xo=np.delete(X,i+1,axis=1); out[c]=1/(1-sm.OLS(yj,Xo).fit().rsquared)
    return out
v = vif(df, ['Knosp','prev_surgery','log2vol'])
coll = pd.DataFrame([{'cohort':'full_n100','n':len(df),
                      'spearman_Knosp_log2vol':round(float(sp.statistic),3),
                      'VIF_Knosp':round(v['Knosp'],2),'VIF_prev_surgery':round(v['prev_surgery'],2),
                      'VIF_log2vol':round(v['log2vol'],2),
                      'interpretation':'Model A=total Knosp effect; Model B=Knosp effect at equal volume (distinct estimands, not a VIF/p-value deletion)'}])
coll.to_csv(RES/'d3_collinearity.csv', index=False)
print('\n=== collinearity (full) ===\n', coll.to_string(index=False))

# ---------- bootstrap coefficient stability (CORRECTED Firth) ----------
def boot(df, y, cols, B=1000, seed=20260925):
    rs = np.random.default_rng(seed)
    Xc = np.column_stack([np.ones(len(df)), df[cols].astype(float).to_numpy()]); yv = df[y].to_numpy().astype(float)
    out=[]
    for _ in range(B):
        idx = rs.integers(0,len(df),len(df))
        try:
            if len(np.unique(yv[idx]))<2: continue
            out.append(firth_fit(Xc[idx], yv[idx]))
        except Exception:
            continue
    out=np.array(out); d={'outcome':y,'cohort':'full_n100','method':'Firth','B_success':int(len(out)),'B_target':B}
    for i,c in enumerate(cols):
        d[f'boot_{c}_median_OR']=round(float(np.exp(np.median(out[:,i+1]))),3)
        d[f'boot_{c}_pct2.5-97.5']=f"{np.exp(np.percentile(out[:,i+1],2.5)):.2f}-{np.exp(np.percentile(out[:,i+1],97.5)):.2f}"
        d[f'boot_{c}_prop_positive']=round(float(np.mean(out[:,i+1]>0)),3)
    return d
bs = pd.DataFrame([boot(frame('residual','full'),'residual',['Knosp','prev_surgery']),
                   boot(frame('subtotal_eor','full'),'subtotal_eor',['Knosp','prev_surgery'])])
bs.to_csv(RES/'d3_bootstrap_coefficient_stability.csv', index=False)
print('\n=== bootstrap coefficient stability (Firth, full) ===\n', bs.to_string(index=False))
print('\n[07 done]')
