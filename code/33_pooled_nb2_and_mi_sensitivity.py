from pathlib import Path; import numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf, json
from scipy import stats as st; import warnings; warnings.filterwarnings("ignore")
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import DATA_DIR as DATA, RESULTS_DIR as RES
SEED=20260925
d3=pd.read_csv(DATA/"d3_clean.csv"); d3["Knosp"]=pd.to_numeric(d3["Knosp"],errors="coerce"); d4=pd.read_csv(DATA/"d4_clean.csv")
tw=pd.DataFrame({"Knosp":d3["Knosp"].astype(float),"prev":pd.to_numeric(d3["prev_surgery"],errors="coerce").astype(float),
                 "los":pd.to_numeric(d3["los_days"],errors="coerce").astype(float),"cohort":0.0})
sp=pd.DataFrame({"Knosp":d4["knosp"].astype(float),"prev":pd.to_numeric(d4["prev_grouped"],errors="coerce").astype(float),
                 "los":pd.to_numeric(d4["los_days"],errors="coerce").astype(float),"cohort":1.0})
def pool(df): 
    df=df.dropna(subset=["Knosp","prev","los"]).copy(); df["Knosp_x"]=df["Knosp"]*df["cohort"]; df["prev_x"]=df["prev"]*df["cohort"]; return df
def nb2(df,cols):
    X=sm.add_constant(df[cols].astype(float)); m=sm.NegativeBinomial(df["los"].astype(float),X,loglike_method="nb2").fit(disp=0)
    conv=bool(m.mle_retvals.get("converged",False))
    return m,float(m.params.iloc[-1]),len(df),conv
def rec(tag,label,m,alpha,n,conv,term="Knosp_x"):
    ci=m.conf_int(); return {"dataset":tag,"model":label,"est":round(float(np.exp(m.params[term])),3),
        "lo":round(float(np.exp(ci.loc[term,0])),3),"hi":round(float(np.exp(ci.loc[term,1])),3),
        "p":round(float(m.pvalues[term]),3),"n":n,"alpha":round(alpha,3),"converged":conv,"tw_gt30_excluded":tag.startswith("EXCL")}
full=pool(pd.concat([tw,sp],ignore_index=True)); ex=pool(pd.concat([tw[tw.los<=30],sp],ignore_index=True))
rows=[]
for tag,df in [("FULL",full),("EXCL_TW_gt30",ex)]:
    mA,aA,nA,cA=nb2(df,["Knosp","prev","cohort","Knosp_x"]); rows.append(rec(tag,"A: +Knosp*cohort",mA,aA,nA,cA))
    mB,aB,nB,cB=nb2(df,["Knosp","prev","cohort","Knosp_x","prev_x"]); rows.append(rec(tag,"B: +prev*cohort too",mB,aB,nB,cB))
    X=sm.add_constant(df[["Knosp","prev","cohort","Knosp_x","prev_x"]].astype(float))
    mP=sm.GLM(df["los"].astype(float),X,family=sm.families.Poisson()).fit(cov_type="HC1"); cP=mP.conf_int()
    rows.append({"dataset":tag,"model":"D: Poisson-HC1 +both","est":round(float(np.exp(mP.params["Knosp_x"])),3),
        "lo":round(float(np.exp(cP.loc["Knosp_x",0])),3),"hi":round(float(np.exp(cP.loc["Knosp_x",1])),3),
        "p":round(float(mP.pvalues["Knosp_x"]),3),"n":len(df),"alpha":None,"converged":bool(mP.converged),"tw_gt30_excluded":tag.startswith("EXCL")})
    # C: separate-cohort NB2, each own alpha; difference Wald z (SE = sqrt(var_tw + var_sp), independent cohorts)
    def kn(d):
        X=sm.add_constant(d[["Knosp","prev"]].astype(float)); m=sm.NegativeBinomial(d["los"].astype(float),X,loglike_method="nb2").fit(disp=0)
        conv=bool(m.mle_retvals["converged"])
        return float(m.params["Knosp"]),float(m.bse["Knosp"]),float(m.params["alpha"]),conv
    twd=tw.dropna(subset=["Knosp","prev","los"])
    if tag.startswith("EXCL"):
        twd=tw[tw.los<=30].dropna(subset=["Knosp","prev","los"])
    spd=sp.dropna(subset=["Knosp","prev","los"])
    b1,s1,a1,c1=kn(twd); b2,s2,a2,c2=kn(spd); z=(b2-b1)/np.sqrt(s1**2+s2**2); pv=2*st.norm.sf(abs(z))
    d=b2-b1; sed=np.sqrt(s1**2+s2**2); rows.append({"dataset":tag,"model":"C: separate NB2 (diff z-test)",
        "est":round(float(np.exp(d)),3),"lo":round(float(np.exp(d-1.96*sed)),3),"hi":round(float(np.exp(d+1.96*sed)),3),
        "p":round(float(pv),3),"n":len(twd)+len(spd),"alpha":f"tw={round(a1,3)},sp={round(a2,3)}",
        "converged":bool(c1 and c2),"tw_gt30_excluded":tag.startswith("EXCL"),"n_tw":len(twd),"n_sp":len(spd)})
df_out=pd.DataFrame(rows)
df_out.to_csv(RES/"sr_pooled_nb2_sensitivity.csv",index=False)
print("=== P0-1 pooled-NB2 interaction sensitivity ==="); print(df_out.to_string(index=False))

# NOTE: An earlier draft included a third (French) cohort and a multiple-imputation
# log(LOS) linearity sensitivity for it (sr_mi_linearity_sensitivity.json). That cohort was
# removed at revision, so this release contains only the two-cohort analyses below.

# P1-4b: de-identified clinical context of the two extreme (>30 d) Taiwanese stays
cols=["Knosp","prev_surgery","los_days","OP time","EOR","Intra-OP CSF","CSF Leakage","ICH","DI (>1m)","2nd OP","other complication","Size"]
ext=d3[pd.to_numeric(d3["los_days"],errors="coerce")>30][cols].copy()
ext.columns=["knosp","prev_surgery","los_days","op_time_min","extent_of_resection","intraop_csf","csf_leakage","intracranial_haemorrhage","diabetes_insipidus_gt1m","second_operation","other_complication","tumour_size"]
ext.to_csv(RES/"sr_los_extreme_clinical_context.csv",index=False)
print("\n=== P1-4b extreme (>30 d) Taiwanese stays, clinical context (de-identified) ==="); print(ext.to_string(index=False))
