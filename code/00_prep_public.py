"""
00_prep_public.py -- France-free data preparation for the two-cohort JNS release.

Builds the Taiwan (index cohort) analytic table from the raw PUBLIC PLOS ONE workbook,
reproducing the cleaning logic of the original 01_prep.py *d3 section only*. The Spain
table is built by 08_d4_prep.py (reads the raw Mendeley .sav). The French cohort (D1)
was removed at revision and is intentionally not part of this release.

Output: data/d3_clean.csv (patient-level; written locally, NOT committed to the repo).

NOTE (portability adaptation vs the original working scripts):
  - Paths are resolved repo-relative via _paths (see data/DATA_MANIFEST.md).
  - The France (d1) branch of the original 01_prep.py was dropped; no statistical
    logic for the two retained cohorts was changed.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import DATA_DIR, RAW_TAIWAN_XLSX

OUT = DATA_DIR
OUT.mkdir(parents=True, exist_ok=True)

# ---------------- D3 (Taiwan, PONE, n=103 raw) ----------------
f = RAW_TAIWAN_XLSX
d3_one = pd.read_excel(f, sheet_name='One hand');  d3_one['approach'] = 'one_hand'
d3_two = pd.read_excel(f, sheet_name='Two Hand');  d3_two['approach'] = 'two_hand'
d3 = pd.concat([d3_one, d3_two], ignore_index=True)

d3 = d3.rename(columns={
    'Previous OP':      'prev_surgery',
    'Tumor_volume':     'volume_cm3',
    'MRI Residual':     'residual',
    'Hospital stay':    'los_days',
    '2 hand':           'two_hand_flag',
})

# EOR categories -> binary outcomes
d3['subtotal_eor'] = d3['EOR'].map({'Subtotal': 1, 'Near Total': 0, 'Total': 0})
d3['incomplete_eor'] = d3['EOR'].map({'Subtotal': 1, 'Near Total': 1, 'Total': 0})

for c in ['Knosp', 'prev_surgery', 'volume_cm3', 'residual', 'subtotal_eor',
          'incomplete_eor', 'los_days', 'Age', 'Hardy']:
    d3[c] = pd.to_numeric(d3[c], errors='coerce')

# Sex / gender: propagate a harmonised 1=male, 0=female column so the age/sex
# sensitivity Firth models (manuscript Supplementary Table S8) are reproducible
# from this pipeline. The PLOS ONE workbook carries sex under one of several
# plausible headers; we take the first that is present.
_SEX_ALIASES = ['Sex', 'sex', 'GENDER', 'Gender', 'gender', 'M/F', 'Gen']
_raw_sex = next((c for c in _SEX_ALIASES if c in d3.columns), None)
if _raw_sex is not None:
    _s = d3[_raw_sex]
    if _s.dtype == object:
        _s = _s.astype(str).str.strip().str.upper()
        d3['male'] = np.where(_s.isin(['M', 'MALE', '1', '1.0']), 1.0,
                     np.where(_s.isin(['F', 'FEMALE', '0', '0.0']), 0.0, np.nan))
    else:
        _sn = pd.to_numeric(_s, errors='coerce')
        d3['male'] = np.where(_sn.isin([1]), 1.0, np.where(_sn.isin([0]), 0.0, np.nan))
    print(f'Sex column detected: {_raw_sex!r}; male n = {int(np.nansum(d3["male"]))}')
else:
    d3['male'] = np.nan
    print('WARNING: no sex/gender column found in raw Taiwan workbook; '
          'd3_clean.csv `male` will be NaN and the D3 age/sex Firth model in '
          'code/44_age_sex_and_threshold_sensitivity.R will be skipped.')

d3['log2vol'] = np.log2(d3['volume_cm3'])
d3['Knosp_high'] = np.where(d3['Knosp'].isna(), np.nan, (d3['Knosp'] >= 3).astype(float))

def knosp_3cat(x):
    if pd.isna(x):
        return np.nan
    if x <= 1:
        return 0.0
    if x == 2:
        return 1.0
    return 2.0

d3['Knosp_3cat'] = d3['Knosp'].map(knosp_3cat)
d3['los_gt7'] = np.where(d3['los_days'].isna(), np.nan, (d3['los_days'] > 7).astype(float))

d3.to_csv(OUT / 'd3_clean.csv', index=False)

# ---------------- sanity summary ----------------
print('=== D3 (Taiwan, PONE) ===')
print('n rows:', len(d3))
print('MRI Residual=1 :', int(d3['residual'].sum()), ' / missing', int(d3['residual'].isna().sum()))
print('Subtotal EOR=1 :', int(d3['subtotal_eor'].sum()))
print('Incomplete=1   :', int(d3['incomplete_eor'].sum()))
print('Knosp dist     :\n', d3['Knosp'].value_counts(dropna=False).sort_index().to_string())
print('prev_surgery=1 :', int(d3['prev_surgery'].sum()))
print('Wrote', str(OUT / 'd3_clean.csv'))
