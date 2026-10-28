"""
09_d4_domain_separation.py — Within-D4 test of the two-domain thesis:
in the SAME 100 patients, does preoperative anatomy predict (a) resectability and
(b) hospitalisation? Expect anatomy->resection strong, anatomy->LOS null.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from _paths import DATA_DIR as _DATA, RESULTS_DIR as _RES
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(_DATA); RES = Path(_RES)
d4 = pd.read_csv(DATA/'d4_clean.csv')

rows = []
# (a) anatomy -> resection (from the Firth R model; recomputed here as MLE for the table)
mA = smf.glm('residual ~ knosp + log2vol', data=d4, family=sm.families.Binomial()).fit()
ci = mA.conf_int()
rows.append({'domain':'Anatomy -> resection (residual)', 'term':'Knosp',
             'est': float(np.exp(mA.params["knosp"])),
             'lo': float(np.exp(ci.loc["knosp",0])), 'hi': float(np.exp(ci.loc["knosp",1])),
             'p': float(mA.pvalues["knosp"]), 'measure':'OR'})
# (b) anatomy -> hospitalisation (continuous LOS, NB2 log-link)
mP = smf.glm('los_days ~ knosp + log2vol', data=d4, family=sm.families.Poisson()).fit()
alpha = max(float(np.sum(mP.resid_pearson**2)/(len(d4)-len(mP.params))), 0.05)
mNB = smf.glm('los_days ~ knosp + log2vol', data=d4, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
ci2 = mNB.conf_int()
for t, lbl in [('knosp','Knosp'), ('log2vol','log2(volume)')]:
    rows.append({'domain':'Anatomy -> hospitalisation (LOS days, NB2)', 'term':lbl,
                 'est': float(np.exp(mNB.params[t])),
                 'lo': float(np.exp(ci2.loc[t,0])), 'hi': float(np.exp(ci2.loc[t,1])),
                 'p': float(mNB.pvalues[t]), 'measure':'IRR'})
# (c) anatomy -> pLOS (>4d) logistic
mL = smf.glm('plos_gt4 ~ knosp + log2vol', data=d4, family=sm.families.Binomial()).fit()
ci3 = mL.conf_int()
rows.append({'domain':'Anatomy -> pLOS (>4d)', 'term':'Knosp',
             'est': float(np.exp(mL.params["knosp"])),
             'lo': float(np.exp(ci3.loc["knosp",0])), 'hi': float(np.exp(ci3.loc["knosp",1])),
             'p': float(mL.pvalues["knosp"]), 'measure':'OR'})
sep = pd.DataFrame(rows)
sep.to_csv(RES/'d4_domain_separation.csv', index=False)
print('=== D4 domain separation ===')
print(sep.round(3).to_string(index=False))

# descriptive
desc = d4.groupby("knosp").agg(n=('los_days','size'), mean_LOS=('los_days','mean'),
                               residual_rate=('residual','mean')).round(2)
desc.to_csv(RES/'d4_los_by_knosp.csv')
print('\n=== mean LOS & residual rate by Knosp (D4) ===')
print(desc.to_string())
by_gtr = d4.groupby('gtr').agg(n=('los_days','size'), mean_LOS=('los_days','mean')).round(2)
print('\n=== mean LOS by GTR (D4) ==='); print(by_gtr.to_string())
by_gtr.to_csv(RES/'d4_los_by_gtr.csv')

# ---- D4 Knosp-only sensitivity (n=100, avoids the 5 ambiguous previous-surgery codes) ----
import sys as _sys; _sys.path.insert(0, str(ROOT/'code'))
from firth_utils import firth_summary
X = np.column_stack([np.ones(len(d4)), d4['knosp'].to_numpy(float)])
tab, _ = firth_summary(X, d4['residual'].to_numpy(float), ['Int','knosp'])
ko = tab[tab.term=='knosp'].iloc[0]
knosp_only = pd.DataFrame([{'model':'D4_Knosp_only_n100','n':len(d4),'term':'knosp',
                            'OR':round(float(ko['OR']),3),'OR_lo':round(float(ko['OR_lo']),3),
                            'OR_hi':round(float(ko['OR_hi']),3),'p':float(ko['p_lr']),
                            'ci_method':'profile penalized-likelihood'}])
knosp_only.to_csv(RES/'d4_knosp_only_sensitivity.csv', index=False)
print('\n=== D4 Knosp-only (n=100) sensitivity ===')
print(knosp_only.to_string(index=False))

# ---- D3 Knosp -> LOS cross-check (D3 raw also carries hospital stay) ----
d3 = pd.read_csv(DATA/'d3_clean.csv').dropna(subset=['Knosp','los_days'])
mP = smf.glm('los_days ~ Knosp + log2vol', data=d3, family=sm.families.Poisson()).fit()
a3 = max(float(np.sum(mP.resid_pearson**2)/(len(d3)-len(mP.params))), 0.05)
mNB3 = smf.glm('los_days ~ Knosp + log2vol', data=d3, family=sm.families.NegativeBinomial(alpha=a3)).fit()
ci3b = mNB3.conf_int()
d3['plos7'] = (d3['los_days']>7).astype(int)
mL7 = smf.glm('plos7 ~ Knosp', data=d3, family=sm.families.Binomial()).fit()
ci7 = mL7.conf_int()
d3los = pd.DataFrame([
 {'analysis':'D3 LOS(NB2) ~ Knosp','term':'Knosp','est':round(float(np.exp(mNB3.params['Knosp'])),3),
  'lo':round(float(np.exp(ci3b.loc['Knosp',0])),3),'hi':round(float(np.exp(ci3b.loc['Knosp',1])),3),'p':round(float(mNB3.pvalues['Knosp']),3),'measure':'IRR'},
 {'analysis':'D3 pLOS>7d ~ Knosp','term':'Knosp','est':round(float(np.exp(mL7.params['Knosp'])),3),
  'lo':round(float(np.exp(ci7.loc['Knosp',0])),3),'hi':round(float(np.exp(ci7.loc['Knosp',1])),3),'p':round(float(mL7.pvalues['Knosp']),3),'measure':'OR'},
])
d3los.to_csv(RES/'d3_knosp_los_crosscheck.csv', index=False)
print('\n=== D3 Knosp -> LOS cross-check (anatomy->LOS is weak/inconsistent, not clean null) ===')
print(d3los.to_string(index=False))
print('mean LOS by Knosp (D3):', d3.groupby('Knosp')['los_days'].mean().round(1).to_dict())
