"""Regenerate France-free, reviewer-simplified Figures 1-3 for JNS v4.
Fig1: cleaner two-cohort / two-lane workflow (shorter labels, larger font).
Fig2: Knosp-grade absolute risk bars (per-grade counts) + adjusted-OR forest; Spain grade-0 n=1 flagged.
Fig3: main panel A = NB2 per-cohort Knosp->LOS only (log-OLS/Poisson moved to Supplementary); panel B = NB2 interaction influence.
Forest numbers read from canonical result CSVs (no hardcoded effect sizes)."""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import ScalarFormatter, FixedLocator
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import RESULTS_DIR, FIG_DIR

RES = str(RESULTS_DIR)
FIG = str(FIG_DIR)
NAVY,BLUE,TEAL,RED,GREY,LGREY = "#3C5488","#4DBBD5","#00A087","#DC0000","#5A5A5A","#D9D9D9"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":9,"axes.linewidth":0.8,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"ps.fonttype":42})
def save(fig,name):
    fig.savefig(f"{FIG}/{name}.png",dpi=300,bbox_inches="tight",facecolor="white"); plt.close(fig); print("wrote",name)

def wilson(k,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=(z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/d; return (c-h,c+h)

def forest(ax,rows,xlabel,xmax,null=1.0,fs=8.5,valpos="right"):
    ys=np.arange(len(rows))[::-1]
    for y,(lab,est,lo,hi,col,mk) in zip(rows and ys,rows):
        ax.plot([lo,hi],[y,y],color=col,lw=1.8,solid_capstyle="round",zorder=2)
        ax.plot(est,y,marker=mk,ms=6,mfc=col,mec=col,zorder=3)
        ax.annotate(f"{est:.2f} ({lo:.2f}\u2013{hi:.2f})",xy=(1.02,y),xycoords=("axes fraction","data"),
                    ha="left",va="center",fontsize=fs-1.5,color=col) if valpos=="right" else \
        ax.text(est,y+0.28,f"{est:.2f} ({lo:.2f}\u2013{hi:.2f})",ha="center",va="bottom",fontsize=fs-1.5,color=col)
    ax.axvline(null,color=GREY,lw=0.8,ls="--",zorder=1); ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows],fontsize=fs)
    for l in ax.get_yticklabels(): l.set_verticalalignment("center")
    ax.set_ylim(-0.7,len(rows)-0.3); ax.set_xscale("log"); ax.set_xlim(0.45,xmax)
    ticks=[t for t in [0.5,1,2,5,10,20,50,100] if 0.45<=t<=xmax]; ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(ScalarFormatter()); ax.set_xticklabels(["%g"%t for t in ticks])
    ax.set_xlabel(xlabel,fontsize=fs); ax.grid(axis="x",color=LGREY,lw=0.5,zorder=0); ax.set_axisbelow(True)

# ---------------- Figure 1 (simplified) ----------------
fig=plt.figure(figsize=(7.0,4.6)); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")
def box(x,y,w,h,t,fc="white",ec=NAVY,fs=9,tc="black",weight="normal",lw=1.2):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.4,rounding_size=1.8",lw=lw,ec=ec,fc=fc))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=fs,color=tc,weight=weight,linespacing=1.3)
def arrow(x1,y1,x2,y2,c=GREY): ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=12,lw=1.4,color=c,shrinkA=0,shrinkB=0))
box(4,88,92,8,"Knosp grade and incomplete resection across two independent pituitary cohorts",fc=NAVY,ec=NAVY,tc="white",fs=10,weight="bold")
box(4,78,92,7,"Secondary analysis of public data (STROBE) — testing reproducibility of association, not validating a prediction model",fc="#EEF3F8",ec=NAVY,fs=8.6)
box(8,56,40,16,"TAIWAN (index)\n103 -> 100 analysed\nKnosp · MRI residual\nLOS · previous surgery",fc="#E7F5F9",ec=BLUE,fs=8.8,weight="bold")
box(52,56,40,16,"SPAIN (replication)\n100 analysed\nKnosp · volumetric GTR\nLOS · previous surgery",fc="#E7F5F9",ec=BLUE,fs=8.8,weight="bold")
box(4,45,92,7,"Shared specification: ordinal Knosp + previous surgery; Firth logistic; harmonisation & leakage check",fc="#FBEEEC",ec=RED,fs=8.6)
arrow(28,55.5,28,52.5); arrow(72,55.5,72,52.5); arrow(50,44.5,50,41.5)
box(8,25,40,15,"PRIMARY — incomplete resection\ndirection reproduced across cohorts\nabsolute risk NOT shown transportable",fc="#EAF6F2",ec=TEAL,fs=8.8,weight="bold")
box(52,25,40,15,"EXPLORATORY — length of stay\nKnosp × cohort interaction\ninfluence of two extreme stays",fc="#EAF6F2",ec=TEAL,fs=8.8,weight="bold")
box(4,12,92,9,"Reproducible association direction ≠ transportable absolute risk ≠ stable LOS association",fc="#EEF3F8",ec=NAVY,fs=9,weight="bold")
box(20,3,20,6,"→ Figure 2",fc="white",ec=GREY,fs=8.6,tc=NAVY,weight="bold"); box(60,3,20,6,"→ Figure 3",fc="white",ec=GREY,fs=8.6,tc=NAVY,weight="bold")
save(fig,"fig1")

# ---------------- Figure 2 ----------------
d3=pd.read_csv(RES+"/firth_d3_logistf.csv"); d4=pd.read_csv(RES+"/firth_d4_logistf.csv"); fp=pd.read_csv(RES+"/firth_pooled_logistf.csv")
tw=[(0,12),(2,24),(1,21),(14,28),(12,15)]; sp=[(0,1),(3,23),(9,31),(17,23),(22,22)]
fig,(a2a,a2b,a2c)=plt.subplots(1,3,figsize=(7.3,2.9),width_ratios=[1,1,1.2],gridspec_kw={"wspace":0.7})
for ax,counts,title,hl in [(a2a,tw,"A  Taiwan · MRI residual",None),(a2b,sp,"B  Spain · incomplete resection",0)]:
    ks=np.arange(5); rates=[k/n*100 for k,n in counts]; ehi=[]; elo=[]
    for k,n in counts:
        lo,hi=wilson(k,n); elo.append((k/n-lo)*100); ehi.append((hi-k/n)*100)
    cols=[BLUE]*5
    if hl==0: cols[0]=LGREY
    ax.bar(ks,rates,width=0.7,color=cols,edgecolor=NAVY,lw=0.7,zorder=3)
    ax.errorbar(ks,rates,yerr=[elo,ehi],fmt="none",ecolor=NAVY,elinewidth=0.9,capsize=2.5,zorder=4)
    for x,(k,n) in enumerate(counts): ax.text(x,min((k/n*100)+wilson(k,n)[1]*100-k/n*100+4,108),f"{k}/{n}",ha="center",va="bottom",fontsize=6.4,color=NAVY)
    ax.set_xticks(ks); ax.set_xticklabels(["0","1","2","3","4"],fontsize=8); ax.set_xlabel("Knosp grade",fontsize=8)
    ax.set_ylabel("Observed rate (%)",fontsize=8); ax.set_ylim(0,116); ax.set_title(title,fontsize=8.4,loc="left",weight="bold",color=NAVY)
    if hl==0: a2b.text(0,min(6,wilson(0,1)[1]*100+4), "n=1", ha="center", va="bottom", fontsize=6.4, color=RED, weight="bold")
rows=[]
for lab,df,model,term,col,mk in [("Taiwan\nModel A",d3,"D3_ModelA_residual","Knosp",TEAL,"o"),
                                  ("Spain\nModel A",d4,"D4_ModelA_grouped_n100","knosp",TEAL,"o"),
                                  ("Pooled\ncommon slope",fp,"POOLED_common_slope","knosp",NAVY,"s"),
                                  ("Interaction\nTW vs SP",fp,"POOLED_interaction","knosp:cohortSpain",RED,"D")]:
    r=df[(df.model==model)&(df.term==term)].iloc[0]; rows.append((lab,r.OR,r.OR_lo,r.OR_hi,col,mk))
forest(a2c,rows,"Adj. OR per Knosp grade (log)",70)
a2c.set_title("C  Independent replication",fontsize=8.4,loc="left",weight="bold",color=NAVY)
fig.text(0.01,-0.09,"Spain grade 0 = single patient (n=1): lowest-grade estimate is imprecise.",fontsize=7,color=RED)
save(fig,"fig2")

# ---------------- Figure 3 (NB2 only + influence) ----------------
los=pd.read_csv(RES+"/continuous_los_robustness_v3.csv"); infl=pd.read_csv(RES+"/sr_los_interaction_sensitivity.csv")
fig,(a3a,a3b)=plt.subplots(1,2,figsize=(7.3,3.0),gridspec_kw={"wspace":0.55,"width_ratios":[1,1]})
r3=[]
for coh,col,mk in [("Taiwan",TEAL,"o"),("Spain",RED,"_")]:
    r=los[(los.cohort==coh)&(los.model=="NB2 joint MLE")&(los.adjustment.str.contains("primary"))].iloc[0]
    r3.append((f"{coh} · NB2",r.est,r.lo,r.hi,col,mk))
forest(a3a,r3,"Knosp effect on continuous LOS (IRR, log)",1.6)
a3a.set_title("A  Cross-cohort LOS (NB2 primary)",fontsize=8.4,loc="left",weight="bold",color=NAVY)
a3a.text(0.02,0.06,"log-OLS HC3 and Poisson-HC1 robustness: Supplementary Table S4",transform=a3a.transAxes,fontsize=6.6,color=GREY)
inb=infl[infl.model=="NB2 interaction"]
order=["FULL (n=200)","Taiwan drop 1 longest (90d)","Taiwan LOS>30d excluded"]
r3b=[( {0:"Full data",1:"Drop 1 longest",2:"Excl. LOS>30d"}[i], inb[inb.dataset==o].iloc[0].est, inb[inb.dataset==o].iloc[0].lo, inb[inb.dataset==o].iloc[0].hi, (TEAL if inb[inb.dataset==o].iloc[0].est<1 else GREY), "o") for i,o in enumerate(order)]
forest(a3b,r3b,"Knosp × cohort ratio of IRRs (log)",1.3)
a3b.set_title("B  Influence analysis",fontsize=8.4,loc="left",weight="bold",color=NAVY)
a3b.text(0.98,0.05,"P: 0.008 → 0.016 → 0.17",transform=a3b.transAxes,ha="right",fontsize=6.8,color=RED)
save(fig,"fig3")
