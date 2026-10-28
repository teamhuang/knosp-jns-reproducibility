"""
CORRECTED adjusted absolute risk (reviewer P0 estimand fix).
Estimand: model-adjusted AVERAGE risk within each observed Knosp band, standardized to
the whole-cohort previous-surgery distribution.
  For each patient i with their OWN Knosp_i, and standardizing previous surgery to the
  cohort prevalence f:
      p_i = f*expit(b0 + b1*Knosp_i + b2*1) + (1-f)*expit(b0 + b1*Knosp_i + b2*0)
  adj_risk(band) = mean_{i in band} p_i       # individual-level prediction (E[expit(Xb)]),
                                              # NOT expit(mean X)
  adjusted RD = adj_risk(high band) - adj_risk(low band)
Bootstrap: 1000 patient-level resamples; refit Firth each replicate, recompute f and the
standardized band risks and RD; report 2.5/97.5 percentile CIs and success count.
Uses the same firth_utils engine (matches canonical R logistf point estimates <0.1%).
"""
import sys, json
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import D3_CLEAN, D4_CLEAN, RESULTS_DIR
CODE = str(Path(__file__).resolve().parent)
from firth_utils import firth_fit

SEED = 20260925
B = 1000
def expit(x): return 1.0/(1.0+np.exp(-x))

def point(df, kc, pc, yc):
    kn = df[kc].to_numpy(float); pv = df[pc].to_numpy(float); yv = df[yc].to_numpy(float)
    X = np.column_stack([np.ones_like(kn), kn, pv])
    b = firth_fit(X, yv)
    f = pv.mean()
    p_i = f*expit(b[0]+b[1]*kn+b[2]) + (1-f)*expit(b[0]+b[1]*kn)   # prev-standardized, own Knosp
    high = kn >= 3
    r_low  = p_i[~high].mean()
    r_high = p_i[high].mean()
    return dict(r_low=r_low, r_high=r_high, rd=r_high-r_low, n=len(df),
                events=int(yv.sum()), nl=int((~high).sum()), nh=int(high.sum()), f=f,
                knosp_OR=float(np.exp(b[1])), prev_OR=float(np.exp(b[2])))

def boot(df, kc, pc, yc):
    rs = np.random.default_rng(SEED); n=len(df)
    rl, rh, rds = [], [], []
    ok=0
    for _ in range(B):
        b=df.iloc[rs.integers(0,n,n)]
        try: r=point(b,kc,pc,yc)
        except Exception: continue
        if not np.isfinite(r['rd']) or r['nl']==0 or r['nh']==0: continue
        ok+=1; rl.append(r['r_low']); rh.append(r['r_high']); rds.append(r['rd'])
    P=lambda a:(float(np.percentile(a,2.5)),float(np.percentile(a,97.5)))
    return ok, P(rl), P(rh), P(rds)

out={}
tw=pd.read_csv(D3_CLEAN).dropna(subset=['Knosp','prev_surgery','volume_cm3','residual'])  # analytic set n=100
sp=pd.read_csv(D4_CLEAN)                                                                   # analytic set n=100
for name,d,kc,pc,yc in [('Taiwan_residual',tw,'Knosp','prev_surgery','residual'),
                        ('Spain_incomplete',sp,'knosp','prev_grouped','residual')]:
    p=point(d,kc,pc,yc); ok,(ll,lh),(hl,hh),(dl,dh)=boot(d,kc,pc,yc)
    out[name]=dict(
        n=p['n'], events=p['events'], n_low_band=p['nl'], n_high_band=p['nh'],
        prev_prevalence_std=round(p['f'],4), knosp_OR=round(p['knosp_OR'],4), prev_OR=round(p['prev_OR'],4),
        adj_risk_low=round(p['r_low'],4),  adj_risk_low_CI=[round(ll,4),round(lh,4)],
        adj_risk_high=round(p['r_high'],4), adj_risk_high_CI=[round(hl,4),round(hh,4)],
        adj_RD=round(p['rd'],4), adj_RD_CI=[round(dl,4),round(dh,4)],
        boot_success=f"{ok}/{B}",
        estimand="model-adjusted average risk within each Knosp band, previous surgery standardized to cohort prevalence; E[expit(Xb)] over individuals (not expit of band-mean)")
print(json.dumps(out,indent=2,ensure_ascii=False))
with open(str(RESULTS_DIR / 'sr_knosp_adjusted_absolute_risk_v2.json'), 'w', encoding='utf-8') as fh:
    fh.write(json.dumps(out, indent=2, ensure_ascii=False))
