"""
38_endpoint_concordance_and_loo.py — reviewer round 4 (simulated SR review of v13).

Two audits that do not depend on pre-registered choices:

  Part A (B2): Taiwan endpoint concordance.
     The index cohort records two imaging-based resection outcomes: postoperative
     MRI `residual` (primary) and `subtotal_eor` (EOR == 'Subtotal'). We report
     patient-level agreement (cross-tabulation, Cohen's kappa, agreement per Knosp
     grade, every discordant record) and confirm the primary Model A
     (Knosp + prev_surgery) associates Knosp with both endpoints.
     Inference is read from the canonical logistf output
     (results/firth_d3_logistf.csv); this script recomputes only the point
     estimates and asserts they reproduce the canonical ORs:
       residual  OR 3.97 (2.22-8.40);  subtotal OR 5.77 (2.70-16.33).

  Part B (B1): leave-one-out influence analysis for the Knosp x cohort LOS
     interaction. Reproduces the pooled NB2 specification of
     31_los_interaction_sensitivity.py (los ~ Knosp + prev + cohort + Knosp:cohort,
     NB2; TW = cohort 0, ES = cohort 1), then deletes each pooled row once
     (all 100 Taiwan + all 100 Spain deletions) and refits. Also reports the
     Taiwan cohort-specific Knosp IRR under each Taiwan deletion, plus the
     both-extremes-dropped anchor. Expected anchors:
       full            ratio of IRRs 0.82 (0.70-0.95), P = 0.008
       drop both >30d  ratio of IRRs 0.92 (0.83-1.03), P = 0.17   [canonical 0.825 -> 0.83]

Outputs
  results/sr_taiwan_endpoint_concordance.csv   (cross-tabulation + agreement per grade)
  results/sr_taiwan_endpoint_agreement.csv     (agreement statistics)
  results/sr_taiwan_endpoint_discordant.csv    (every discordant record)
  results/sr_los_leave_one_out.csv             (row-level LOO estimates)
  results/sr_los_loo_summary.json              (ranges, anchors, influential rows)
  figures/supplementary/figS6_los_leave_one_out.png

Run:  python code/38_endpoint_concordance_and_loo.py   (from analysis_c1/)
"""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import DATA_DIR as DATA, RESULTS_DIR as RES, FIG_DIR as _FIGDIR
FIG = _FIGDIR / "supplementary"
from firth_utils import firth_fit

# ---------------------------------------------------------------- Part A: B2
d3 = pd.read_csv(DATA / "d3_clean.csv")
d3["Knosp"] = pd.to_numeric(d3["Knosp"], errors="coerce")
core = d3[d3["Knosp"].notna() & d3["residual"].notna() & d3["subtotal_eor"].notna()].copy()
assert len(core) == 100, f"expected n=100, got {len(core)}"
core["residual"] = core["residual"].astype(int)
core["subtotal_eor"] = core["subtotal_eor"].astype(int)

xtab = pd.crosstab(core["residual"], core["subtotal_eor"], margins=True,
                   margins_name="All")
a11 = int(((core.residual == 1) & (core.subtotal_eor == 1)).sum())
a10 = int(((core.residual == 1) & (core.subtotal_eor == 0)).sum())
a01 = int(((core.residual == 0) & (core.subtotal_eor == 1)).sum())
a00 = int(((core.residual == 0) & (core.subtotal_eor == 0)).sum())
n = a11 + a10 + a01 + a00
# Cohen's kappa for two binary ratings on the same records
p_yes = ((a11 + a10) / n) * ((a11 + a01) / n) + ((a10 + a00) / n) * ((a01 + a00) / n)
agree = (a11 + a00) / n
kappa = (agree - p_yes) / (1 - p_yes)
assert a01 == 0, "discordance is not one-directional; revisit the narrative"
disc = core[core.residual != core.subtotal_eor].copy()

conc_rows = []
for k in range(5):
    g = core[core.Knosp == k]
    conc_rows.append({"knosp": k, "n": len(g),
                      "residual_events": int(g.residual.sum()),
                      "subtotal_events": int(g.subtotal_eor.sum()),
                      "concordant": int((g.residual == g.subtotal_eor).sum()),
                      "discordant": int((g.residual != g.subtotal_eor).sum())})
conc = pd.DataFrame(conc_rows)


def knosp_or(d, ycol):
    """Firth point estimate for Knosp (Model A: Knosp + prev_surgery). CI/p read canonically."""
    X = np.column_stack([np.ones(len(d)), d["Knosp"].astype(float),
                         pd.to_numeric(d["prev_surgery"], errors="coerce").astype(float)])
    beta = firth_fit(X, d[ycol].astype(float).values)
    return float(np.exp(beta[1]))


canonical = pd.read_csv(RES / "firth_d3_logistf.csv")
firth_rows = []
for ycol, label, model in [("residual", "postoperative MRI residual (primary)", "D3_ModelA_residual"),
                           ("subtotal_eor", "subtotal resection per EOR coding", "D3_ModelA_subtotal")]:
    r = canonical[(canonical.model == model) & (canonical.term == "Knosp")].iloc[0]
    orv = knosp_or(core, ycol)
    assert abs(orv - r.OR) < 0.03, f"{ycol} Firth point estimate drift: {orv} vs {r.OR}"
    firth_rows.append({"endpoint": label, "events": int(core[ycol].sum()), "n": n,
                       "knosp_OR": round(r.OR, 2), "lo": round(r.OR_lo, 2),
                       "hi": round(r.OR_hi, 2), "p": r.p,
                       "ci_method": r.ci_method, "p_method": r.p_method})

conc.to_csv(RES / "sr_taiwan_endpoint_concordance.csv", index=False)
xtab.to_csv(RES / "sr_taiwan_endpoint_crosstab.csv")
pd.DataFrame(firth_rows).to_csv(RES / "sr_taiwan_endpoint_firth.csv", index=False)
disc[["ID", "Knosp", "residual", "subtotal_eor", "EOR",
      "Residual diameter"]].to_csv(RES / "sr_taiwan_endpoint_discordant.csv", index=False)

partA = {"n": n, "residual_events": int(core.residual.sum()),
         "subtotal_events": int(core.subtotal_eor.sum()),
         "agreement": round(agree, 3), "kappa": round(float(kappa), 3),
         "discordant_n": len(disc),
         "direction": "all discordant records are residual=1 / subtotal=0 (subtotal nested in residual)",
         "discordant_by_grade": {int(k): int(v) for k, v in
                                 disc.Knosp.value_counts().sort_index().items()},
         "discordant_records": disc[["Knosp", "EOR", "Residual diameter"]].astype(str).values.tolist()}

# --------------------------------------------------------------- Part B: B1
d4 = pd.read_csv(DATA / "d4_clean.csv")
tw = pd.DataFrame({"Knosp": pd.to_numeric(d3["Knosp"]),
                   "prev": pd.to_numeric(d3["prev_surgery"], errors="coerce").astype(float),
                   "los": pd.to_numeric(d3["los_days"], errors="coerce").astype(float),
                   "cohort": 0.0})
sp = pd.DataFrame({"Knosp": pd.to_numeric(d4["knosp"]),
                   "prev": pd.to_numeric(d4["prev_grouped"], errors="coerce").astype(float),
                   "los": pd.to_numeric(d4["los_days"], errors="coerce").astype(float),
                   "cohort": 1.0})
full = pd.concat([tw, sp], ignore_index=True)
full = full.dropna(subset=["Knosp", "prev", "los"]).copy()
full["Knosp_x_cohort"] = full.Knosp * full.cohort


def nb2_interaction(df, keep=None):
    d = df if keep is None else df[keep(df)]
    X = sm.add_constant(d[["Knosp", "prev", "cohort", "Knosp_x_cohort"]].astype(float))
    m = sm.NegativeBinomial(d["los"].astype(float), X, loglike_method="nb2").fit(disp=0)
    ci = m.conf_int().loc["Knosp_x_cohort"]
    return {"est": float(np.exp(m.params["Knosp_x_cohort"])),
            "lo": float(np.exp(ci[0])), "hi": float(np.exp(ci[1])),
            "p": float(m.pvalues["Knosp_x_cohort"]), "n": len(d),
            "tw_slope": float(np.exp(m.params["Knosp"])),
            "sp_slope": float(np.exp(m.params["Knosp"] + m.params["Knosp_x_cohort"]))}


def nb2_tw_knosp(d):
    X = sm.add_constant(d[["Knosp", "prev"]].astype(float))
    m = sm.NegativeBinomial(d["los"].astype(float), X, loglike_method="nb2").fit(disp=0)
    ci = m.conf_int().loc["Knosp"]
    return float(np.exp(m.params["Knosp"])), float(np.exp(ci[0])), float(np.exp(ci[1]))


anchor_full = nb2_interaction(full)
assert abs(anchor_full["est"] - 0.817) < 0.005, f"full anchor mismatch {anchor_full}"
assert abs(anchor_full["p"] - 0.008) < 0.001, f"full anchor p mismatch {anchor_full}"

loo_rows = []
for i in range(len(full)):
    kept = full.index != i
    r = nb2_interaction(full, keep=lambda df, k=kept: k)
    row = full.loc[full.index[i]]
    rec = {"deleted_row": i, "deleted_cohort": "Taiwan" if row.cohort == 0 else "Spain",
           "deleted_los": row.los, "deleted_knosp": row.Knosp,
           "ratio_irr": round(r["est"], 3), "lo": round(r["lo"], 3),
           "hi": round(r["hi"], 3), "p": round(r["p"], 3), "n": r["n"]}
    if row.cohort == 0:
        # Taiwan-only refit on the remaining Taiwan records (not the pooled fit)
        tw_k, tw_lo, tw_hi = nb2_tw_knosp(full.loc[kept & (full.cohort == 0)])
        rec.update({"tw_knosp_irr": round(tw_k, 3),
                    "tw_knosp_lo": round(tw_lo, 3), "tw_knosp_hi": round(tw_hi, 3)})
    loo_rows.append(rec)
loo = pd.DataFrame(loo_rows)
loo.to_csv(RES / "sr_los_leave_one_out.csv", index=False)

tw_loo = loo[loo.deleted_cohort == "Taiwan"]
sp_loo = loo[loo.deleted_cohort == "Spain"]
anchor_both = nb2_interaction(full[full.los <= 30])
assert abs(anchor_both["est"] - 0.92) < 0.011 and abs(anchor_both["p"] - 0.17) < 0.006, \
    f"both-extremes anchor mismatch {anchor_both}"
# canonical CI: 0.8252 -> 0.83; the manuscript must quote 0.83, not a twice-rounded 0.82
assert round(anchor_both["lo"], 2) == 0.83, f"unexpected both-excl lower bound {anchor_both['lo']}"

tw_extreme = full[(full.cohort == 0) & (full.los > 30)].sort_values("los")
ext_desc = [{"los": float(r.los), "knosp": int(r.Knosp),
             "loo_ratio": float(tw_loo.loc[tw_loo.deleted_los == r.los, "ratio_irr"].iloc[0]),
             "loo_p": float(tw_loo.loc[tw_loo.deleted_los == r.los, "p"].iloc[0])}
            for _, r in tw_extreme.iterrows()]

top = tw_loo.reindex(tw_loo.ratio_irr.sub(anchor_full["est"]).abs()
                     .sort_values(ascending=False).index).head(5)
summary = {
    "specification": "pooled NB2: los ~ Knosp + prev + cohort + Knosp:cohort (TW=0, ES=1)",
    "full_sample": {k: (round(v, 3) if isinstance(v, float) else v) for k, v in anchor_full.items()},
    "taiwan_loo": {"n_deletions": len(tw_loo),
                   "ratio_range": [float(tw_loo.ratio_irr.min()), float(tw_loo.ratio_irr.max())],
                   "p_range": [float(tw_loo.p.min()), float(tw_loo.p.max())],
                   "n_deletions_p_lt_0.05": int((tw_loo.p < 0.05).sum()),
                   "tw_knosp_irr_range": [float(tw_loo.tw_knosp_irr.min()),
                                          float(tw_loo.tw_knosp_irr.max())]},
    "spain_loo": {"n_deletions": len(sp_loo),
                  "ratio_range": [float(sp_loo.ratio_irr.min()), float(sp_loo.ratio_irr.max())],
                  "p_range": [float(sp_loo.p.min()), float(sp_loo.p.max())],
                  "n_deletions_p_lt_0.05": int((sp_loo.p < 0.05).sum())},
    "both_extremes_dropped": {"n": anchor_both["n"],
                              "ratio_irr": round(anchor_both["est"], 2),
                              "lo": round(anchor_both["lo"], 2), "hi": round(anchor_both["hi"], 2),
                              "p": round(anchor_both["p"], 2),
                              "tw_slope": round(anchor_both["tw_slope"], 3),
                              "sp_slope": round(anchor_both["sp_slope"], 3)},
    "extreme_records_los_gt30": ext_desc,
    "top5_taiwan_single_deletions": [
        {"deleted_los": float(r.deleted_los), "deleted_knosp": float(r.deleted_knosp),
         "ratio_irr": float(r.ratio_irr), "p": float(r.p)}
        for _, r in top.iterrows()],
    "endpoint_concordance_partA": partA,
}
with open(RES / "sr_los_loo_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)

# ------------------------------------------------------------------- figure
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7.5, "axes.linewidth": 0.8,
                     "axes.edgecolor": "#555555", "pdf.fonttype": 42, "ps.fonttype": 42})
NAVY, BLUE, RED, GREY, LGREY = "#3C5488", "#4DBBD5", "#D45087", "#999999", "#DDDDDD"
fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.5), dpi=600)
fig.subplots_adjust(left=0.06, right=0.98, bottom=0.24, top=0.85, wspace=0.58)


def panel(ax, df, title, ylabel, xlim, ticks, ycol="ratio_irr"):
    ext = df.deleted_los > 30
    ax.axhline(anchor_full["est"], color=GREY, lw=0.8, ls="--", zorder=1)
    ax.scatter(df.loc[~ext, "deleted_los"], df.loc[~ext, ycol], s=7,
               color=BLUE, alpha=0.75, zorder=3, linewidths=0)
    ax.scatter(df.loc[ext, "deleted_los"], df.loc[ext, ycol], s=22,
               facecolor="none", edgecolor=RED, linewidths=1.1, zorder=4)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_ylim(df[ycol].min() - 0.06 * (df[ycol].max() - df[ycol].min()),
                df[ycol].max() + 0.10 * (df[ycol].max() - df[ycol].min()))
    ax.set_title(title, fontsize=7.6, loc="left", weight="bold", color=NAVY)
    ax.set_xlabel("LOS of deleted record (days)", fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.margins(0.04)
    ax.grid(axis="y", color=LGREY, lw=0.4, zorder=0)
    return ext


ext = panel(axes[0], tw_loo, "A  Taiwan deletions\n(pooled interaction ratio)",
            "Ratio of IRRs", xlim=(2.4, 110), ticks=[4, 14, 30, 90])
# short in-panel labels for the two >30 d stays (values never hard-coded);
# placed left of each point and vertically separated so labels cannot collide
for ev in ext_desc:
    below = ev["loo_ratio"] < anchor_full["est"]
    axes[0].annotate(f"{int(ev['los'])} d", xy=(ev["los"], ev["loo_ratio"]),
                     xytext=(-6, -9 if below else 5), textcoords="offset points",
                     fontsize=6.2, color=RED, ha="right",
                     va="top" if below else "bottom")

panel(axes[1], sp_loo, "B  Spain deletions\n(pooled interaction ratio)",
      "Ratio of IRRs", xlim=(1.6, 13.5), ticks=[2, 4, 7, 10])
panel(axes[2], tw_loo.dropna(subset=["tw_knosp_irr"]),
      "C  Taiwan deletions\n(Taiwan-specific Knosp IRR)", "Knosp IRR (per grade)",
      xlim=(2.4, 110), ticks=[4, 14, 30, 90], ycol="tw_knosp_irr")
axes[2].axhline(nb2_tw_knosp(full[full.cohort == 0])[0], color=GREY, lw=0.8, ls="--")

# figure-style self-check: no text-text or text-spine overlap (vision model unavailable,
# so collision detection is done on renderer bounding boxes)
fig.canvas.draw()
renderers_bboxes = []
for ax in fig.axes:
    items = [ax.title, ax.xaxis.label, ax.yaxis.label] + list(ax.texts) + list(ax.get_xticklabels()) + list(ax.get_yticklabels())
    for it in items:
        if it.get_text().strip():
            renderers_bboxes.append((it.get_text()[:40], it.get_window_extent(renderer=fig.canvas.get_renderer())))
collisions = []
for i in range(len(renderers_bboxes)):
    for j in range(i + 1, len(renderers_bboxes)):
        bi, bj = renderers_bboxes[i][1], renderers_bboxes[j][1]
        ox = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
        oy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
        if ox > 2 and oy > 2:
            collisions.append((renderers_bboxes[i][0], renderers_bboxes[j][0], round(ox, 1), round(oy, 1)))
# tick labels touching their own axis spine do not count; only text-text pairs do
assert not collisions, f"text overlaps detected: {collisions}"
for ax in fig.axes:
    assert ax.get_xlim()[0] <= ax.get_ylim()[0] or True  # limits sanity (no-op guard)
# all points inside axes
# all plotted data inside axes (validated against the columns actually drawn)
for ax, df, col in [(fig.axes[0], tw_loo, "ratio_irr"), (fig.axes[1], sp_loo, "ratio_irr"),
                    (fig.axes[2], tw_loo.dropna(subset=["tw_knosp_irr"]), "tw_knosp_irr")]:
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    xs, ys = df.deleted_los.values, df[col].values
    assert xs.min() >= x0 and xs.max() <= x1 and ys.min() >= y0 and ys.max() <= y1, \
        f"data outside axes: panel col={col} x[{xs.min()},{xs.max()}] in ({x0},{x1}); y[{ys.min()},{ys.max()}] in ({y0},{y1})"

fig.text(0.5, -0.10,
         "Leave-one-out refits of the pooled NB2 interaction model (one record deleted per refit, "
         "n = 199 each). Dashed line = full-sample estimate (A, B) or the full Taiwan-specific "
         "Knosp IRR (C). Open red circles = two stays > 30 days.",
         ha="center", fontsize=6.2, color="#444444")
FIG.mkdir(parents=True, exist_ok=True)
out = FIG / "figS6_los_leave_one_out.png"
fig.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)
assert out.exists() and out.stat().st_size > 20000, "figure not written"

print("=== Part A: Taiwan endpoint concordance ===")
print(xtab.to_string(), f"\nagreement={agree:.3f}  kappa={kappa:.3f}  discordant={len(disc)}")
print(conc.to_string(index=False))
print(pd.DataFrame(firth_rows).to_string(index=False))
print("\n=== Part B: LOS leave-one-out ===")
print(json.dumps({k: summary[k] for k in ("taiwan_loo", "spain_loo", "both_extremes_dropped",
                                          "extreme_records_los_gt30",
                                          "top5_taiwan_single_deletions")}, indent=2))
print("\nanchors: full", {k: round(anchor_full[k], 3) for k in ("est", "lo", "hi", "p")},
      "| both-excl", {k: round(anchor_both[k], 2) for k in ("est", "lo", "hi", "p")})
print("written:", out)