"""
France-free Figure 1 (two-cohort study-design / analytical workflow schematic).
Matches the style of 34_submission_figures_v11.py (Arial, Nature palette, spine-thin).
Output -> workspace/figs/fig1.png (300+ dpi) for embedding in the JNS v3 draft.
No numeric results are plotted (pure schematic); all box counts (103->100, 100) come from
the frozen provenance manifest.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import FIG_DIR

NAVY, BLUE, TEAL, RED = "#3C5488", "#4DBBD5", "#00A087", "#DC0000"
GREY, LGREY = "#5A5A5A", "#D9D9D9"
plt.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.5, "pdf.fonttype": 42, "ps.fonttype": 42})

fig = plt.figure(figsize=(7.2, 7.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, fc="white", ec=NAVY, lw=1.1, fs=8, tc="black", weight="normal", rad=1.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.4,rounding_size={rad}",
                                linewidth=lw, edgecolor=ec, facecolor=fc, mutation_aspect=1))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, color=tc,
            weight=weight, linespacing=1.35)

def arrow(x1, y1, x2, y2, color=GREY, lw=1.3):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                                 linewidth=lw, color=color, shrinkA=0, shrinkB=0))

box(6, 91, 88, 6.4, "Reproducibility boundaries of Knosp associations across outcome and measurement systems: study design",
    fc=NAVY, ec=NAVY, tc="white", fs=8.4, weight="bold")
box(6, 83, 88, 6.4, "Retrospective secondary analysis · two independent public individual-level pituitary-surgery cohorts\n"
                    "Association-focused (STROBE) — not development or external validation of a prediction model",
    fc="#EEF3F8", ec=NAVY, fs=7.6)

col_y, col_h = 65, 15
box(8, col_y, 39, col_h, "TAIWAN  (index cohort)\n103 raw → 100 analysed\nKnosp · MRI residual\nLOS · previous surgery",
    fc="#E7F5F9", ec=BLUE, fs=7.6, weight="bold")
box(53, col_y, 39, col_h, "SPAIN  (replication cohort)\n100 analysed\nKnosp · volumetric GTR\nLOS · previous surgery",
    fc="#E7F5F9", ec=BLUE, fs=7.6, weight="bold")
for x in (27.5, 72.5):
    arrow(x, col_y - 0.2, x, 60.2)

box(6, 54, 88, 5.6, "Data harmonisation & quality audit: outcome definitions · covariate coding · missingness · temporal-information / leakage check",
    fc="#FBEEEC", ec=RED, fs=7.4)
arrow(50, 53.8, 50, 50.2)

lane_y, lane_h = 33, 15.5
box(8, lane_y, 39, lane_h,
    "PRIMARY analysis\nKnosp → incomplete resection\nTW + SP (n = 200 combined)\n\nCohort-specific Firth logistic +\nprevious surgery;\nindependent replication of direction",
    fc="#EAF6F2", ec=TEAL, fs=7.2, weight="bold")
box(53, lane_y, 39, lane_h,
    "EXPLORATORY analysis\nKnosp → continuous LOS\nTW + SP (n = 200 combined)\n\nNB2 · Knosp × cohort interaction ·\nextreme-stay influence diagnostics",
    fc="#EAF6F2", ec=TEAL, fs=7.2, weight="bold")
for x in (27.5, 72.5):
    arrow(x, 53.8, x, lane_y + lane_h + 0.2)

box(6, 24, 88, 6.0, "Interpretation of reproducibility and its limits", fc=NAVY, ec=NAVY, tc="white", fs=8, weight="bold")
box(6, 13.5, 88, 9.0,
    "Resection association directionally replicated across cohorts despite differing imaging endpoints\n"
    "Model-adjusted absolute risk reported within cohort; transportability of absolute risk NOT claimed\n"
    "LOS cross-cohort heterogeneity influence-sensitive (depends on two extreme complication-associated stays)\n"
    "No prediction-model validation · no causal interpretation",
    fc="#EEF3F8", ec=NAVY, fs=6.9)

box(18, 6.5, 26, 5.5, "→ Figure 2", fc="white", ec=GREY, fs=7.8, weight="bold", tc=NAVY)
box(56, 6.5, 26, 5.5, "→ Figure 3", fc="white", ec=GREY, fs=7.8, weight="bold", tc=NAVY)

fig.savefig(str(FIG_DIR / "fig1.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(str(FIG_DIR / "fig1.pdf"), bbox_inches="tight", facecolor="white")
print("wrote France-free fig1.png/.pdf")
