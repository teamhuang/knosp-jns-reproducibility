"""
08_d4_prep.py — D4 independent-replication cohort (Mendeley Data, DOI 10.17632/kb55bvtv3y.1,
CC BY 4.0; Pérez-López et al., World Neurosurgery 2021, PMID 33309892). 100 non-functioning
pituitary adenoma operations at a SINGLE centre (Hospital Universitario La Paz, Madrid) by
one surgical team. The dataset's `Hospital` field lists 15 named institutions = REFERRING
hospitals, not operating centres. Used to replicate the anatomical->resectability association
and to explore anatomy->LOS.

Cleaning rules (documented, no silent recoding):
  - residual = 1 if GTR == 'Parcial' (incomplete), 0 if 'Completa' (GTR).
  - previous pituitary surgery: use the dataset's own derived binary `IQPrevia_Grupo`
    (0=No,1=Yes), which collapses raw codes 1/2/3 into "yes" (20 yes / 80 no, complete for
    all 100). Codes 2/3 are prior-surgery subcategories per the authors' recode, NOT errors.
  - Knosp 0-4 ordinal (grade 0 n=1 -> sparse; Firth used).
  - volume (mL) -> log2; LOS (days) -> pLOS at >4d and >median.
  - los_days = `QxEstanciaDias` ("Estancia"); prefix Qx = post-Quirúrgico, and ICU stay is a
    separate nested variable, so Estancia is the postoperative hospital length of stay
    (definition inferred from dataset structure; no README provided -> confirm with authors).
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pyreadstat, pandas as pd, numpy as np
from _paths import DATA_DIR as DATA, RESULTS_DIR as RES, RAW_SPAIN_SAV
DATA = Path(DATA); RES = Path(RES)
df, meta = pyreadstat.read_sav(str(RAW_SPAIN_SAV), apply_value_formats=False)

def num(x): return pd.to_numeric(x, errors='coerce')

d4 = pd.DataFrame({
    'hospital':      num(df['Hospital']),                          # referring institution
    'age':           num(df['Edad']),
    'male':          (num(df['Sexo']) == 1).astype(int),          # 1=Varón
    'knosp':         num(df['RxKnosPre']),
    'volume_ml':     num(df['RxVolPre_ml']),
    'residual_ml':   num(df['RxVolPost_ml']),
    'prev_surgery':  num(df['IQPrevia_Grupo']),                   # authors' binary (0/1, complete)
    'gtr':           (num(df['GTR']) == 1).astype(int),           # 1=Completa
    'resection_pct': num(df['Grado_Reseccion']),
    'suprasellar':   (num(df['RxSupraselar']) == 1).astype(int),
    'clivus':        (num(df['RxClivus']) == 1).astype(int),
    'third_ventr':   (num(df['RxIII_Ventr']) == 1).astype(int),
    'apoplexy':      (num(df['Apoplejia']) == 1).astype(int),
    't2_iso':        (num(df['T2_señal']) == 0).astype(int),      # 0=Isointense
    'expanded_sag':  (num(df['QxExpandidoSagital']) == 1).astype(int),
    'expanded_cor':  (num(df['QxExpandidoCoronal']) == 1).astype(int),
    'di':            (num(df['DI']) == 1).astype(int),
    'reintervention':(num(df['SegReintervencion']) == 1).astype(int),
    'icu_days':      num(df['ComplQxDiasUCI']),
    'los_days':      num(df['QxEstanciaDias']),
})
d4['residual'] = 1 - d4['gtr']                                    # 1 = incomplete resection
d4['log2vol'] = np.log2(d4['volume_ml'].clip(lower=0.1))
d4['knosp_high'] = (d4['knosp'] >= 3).astype(int)                 # 3-4 vs 0-2
d4['plos_gt4'] = (d4['los_days'] > 4).astype(int)
d4['plos_gtmedian'] = (d4['los_days'] > d4['los_days'].median()).astype(int)

# Previous-surgery codings (three variants used by the primary model + sensitivity analyses).
# Source column IQ_hipofisis_previa uses codes {0=none, 1=pituitary, 2/3=ambiguous other}.
# IQPrevia_Grupo is the authors' pre-grouped 0/1 (codes 1/2/3 -> yes).
hip = num(df['IQ_hipofisis_previa'])
d4['prev_grouped'] = d4['prev_surgery']                              # PRIMARY: 0/1, 20 yes
d4['prev_strict']  = np.where(hip.isin([2, 3]), np.nan, (hip == 1).astype(float))  # documented only: 95 (15 yes, 5 dropped)
d4['prev_ambl2no'] = (hip == 1).astype(int)                          # ambiguous 2/3 -> No: 100 (15 yes)

d4.to_csv(DATA/'d4_clean.csv', index=False)

# data dictionary
dd = pd.DataFrame([
    ('residual', 'incomplete resection (GTR Parcial)', f"{int(d4.residual.sum())}/{len(d4)}"),
    ('gtr', 'gross-total resection (GTR Completa)', f"{int(d4.gtr.sum())}/{len(d4)}"),
    ('knosp', 'preop Knosp 0-4', d4.knosp.value_counts().sort_index().to_dict().__str__()),
    ('volume_ml', 'preop tumour volume (mL)', f"median {d4.volume_ml.median():.1f}"),
    ('prev_surgery', 'previous pituitary surgery (authors binary IQPrevia_Grupo; raw 1/2/3 -> yes)', f"usable {int(d4.prev_surgery.notna().sum())}, yes {int(d4.prev_surgery.sum())}"),
    ('los_days', 'postop hospital LOS (Estancia; Qx=post-surgical; ICU separate) - definition to confirm', f"median {d4.los_days.median():.0f}, range {d4.los_days.min():.0f}-{d4.los_days.max():.0f}"),
    ('plos_gt4', 'pLOS > 4 days', f"{int(d4.plos_gt4.sum())}/{len(d4)}"),
    ('hospital', 'referring institutions (NOT operating centres; single-centre La Paz)', f"{int(d4.hospital.nunique())} named institutions"),
])
dd.columns = ['variable', 'definition', 'value']
dd.to_csv(RES/'d4_data_dictionary.csv', index=False)
print('=== D4 prepped: n=%d, %d referring institutions (single operating centre) ===' % (len(d4), d4.hospital.nunique()))
print(dd.to_string(index=False))
