"""Repo-relative path resolution for the Knosp-JNS two-cohort release.

Every analysis script imports this so the pipeline runs from any checkout without
editing hardcoded paths. All directories can be overridden by environment variables:

    KNOSP_DATA_DIR       -> directory holding the two raw PUBLIC files (default: <repo>/data)
    KNOSP_RESULTS_DIR    -> where model-output tables are written      (default: <repo>/results)
    KNOSP_FIG_DIR        -> where figures are written                  (default: <repo>/figures)
    KNOSP_TAIWAN_XLSX    -> explicit path to the Taiwan raw workbook
    KNOSP_SPAIN_SAV      -> explicit path to the Spain raw .sav

Raw PUBLIC input files expected in DATA_DIR (see data/DATA_MANIFEST.md):
    PONE-D-21-14669R1-Raw_data.xlsx   (Taiwan / index cohort, PLOS ONE figshare supplement)
    D4_mendeley_nfpa.sav              (Spain  / replication cohort, Mendeley Data NFPA)
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

DATA_DIR = Path(os.environ.get("KNOSP_DATA_DIR", str(REPO / "data")))
RESULTS_DIR = Path(os.environ.get("KNOSP_RESULTS_DIR", str(REPO / "results")))
FIG_DIR = Path(os.environ.get("KNOSP_FIG_DIR", str(REPO / "figures")))

RAW_TAIWAN_XLSX = Path(os.environ.get(
    "KNOSP_TAIWAN_XLSX", str(DATA_DIR / "PONE-D-21-14669R1-Raw_data.xlsx")))
RAW_SPAIN_SAV = Path(os.environ.get(
    "KNOSP_SPAIN_SAV", str(DATA_DIR / "D4_mendeley_nfpa.sav")))

# Derived (patient-level) analytic tables -- regenerated locally from raw, NOT redistributed.
D3_CLEAN = DATA_DIR / "d3_clean.csv"
D4_CLEAN = DATA_DIR / "d4_clean.csv"

for _d in (RESULTS_DIR, FIG_DIR, FIG_DIR / "supplementary"):
    _d.mkdir(parents=True, exist_ok=True)
