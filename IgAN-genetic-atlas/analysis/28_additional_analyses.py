"""Step 28. Additional analyses, computed from Supplementary Data S1 (Tables S1-S8).
Usage: python 28_additional_analyses.py Supplementary_Data_S1_v1.1.3.xlsx additional_analyses.xlsx
Also writes Aim1_decomposition_bootstrapCI.tsv (used by figures/21_fig2_ancestry.R).

1. Prior sensitivity of coloc.abf PP.H4 to p12 (exact: only H4 depends on p12).
2. Immune vs kidney window counts under different p12.
3. CFH-excluded versions of the primary count, matched-testing analyses
   (exact enumeration instead of Monte Carlo) and the gene-level comparison.
4. Parametric-bootstrap 95% CIs for the ancestry variance decomposition.
5. Steiger directionality for the single-instrument Wald ratios.
"""
import sys
from itertools import combinations
from math import comb, log10
import numpy as np
import pandas as pd
from scipy import stats

X = pd.read_excel(sys.argv[1], sheet_name=None)
OUT = sys.argv[2]
S1, S2, S4, S5, S6 = X["S1"], X["S2"], X["S4"], X["S5"], X["S6"]
strata = X["S8_strata"]
out = {}


def rescale(pp4, r):
    pp4 = np.asarray(pp4, float)
    return r * pp4 / (r * pp4 + (1 - pp4))


P12 = [1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4]

# ---------- 1. prior sensitivity for the main findings ----------
def best(df, **kw):
    d = df
    for k, v in kw.items():
        d = d[d[k] == v]
    return d.sort_values("PP.H4", ascending=False).iloc[0]

rows = []
main = [
    ("TNFSF13", "Classical monocytes (ImmuNexUT)", "EAS", best(S1, gene="TNFSF13", cell="CL_Mono", ancestry="EAS")),
    ("TNFSF13", "Myeloid DC (ImmuNexUT)", "EAS", best(S1, gene="TNFSF13", cell="mDC", ancestry="EAS")),
    ("TNFSF13", "Whole-blood mRNA (JCTF)", "EAS", best(S5, gene="TNFSF13", layer="eQTL", ancestry="EAS")),
    ("APRIL", "Plasma protein (JCTF)", "EAS", best(S5, gene="TNFSF13", layer="pQTL", ancestry="EAS")),
    ("FCAR", "Junction 54,885,525-54,885,805, monocytes (Quach)", "EAS", best(S4, gene="FCAR", study="Quach_2016", quant="leafcutter", ancestry="EAS")),
    ("FCAR", "Junction 54,885,525-54,885,805, whole blood (GTEx)", "EAS", best(S4, gene="FCAR", study="GTEx", tissue="blood", quant="leafcutter", ancestry="EAS")),
    ("Factor H", "Plasma protein (JCTF)", "EUR", best(S5, gene="CFH", layer="pQTL", ancestry="EUR")),
    ("Factor H", "Plasma protein (JCTF)", "EAS", best(S5, gene="CFH", layer="pQTL", ancestry="EAS")),
    ("FHR-5", "Plasma protein (JCTF)", "EUR", best(S5, gene="CFHR5", layer="pQTL", ancestry="EUR")),
    ("CFHR3", "Liver (GTEx)", "EAS", best(S4, gene="CFHR3", tissue="liver", ancestry="EAS")),
    ("CFHR4", "Liver (GTEx)", "EUR", best(S4, gene="CFHR4", tissue="liver", ancestry="EUR")),
    ("CARD9", "Monocytes (ImmuNexUT)", "EUR", best(S1, gene="CARD9", ancestry="EUR")),
    ("ITGAM", "Classical monocytes (ImmuNexUT)", "EAS", best(S1, gene="ITGAM", cell="CL_Mono", ancestry="EAS")),
    ("CD28", "Naive CD4+ T (ImmuNexUT)", "EAS", best(S1, gene="CD28", ancestry="EAS")),
    ("LYN", "B-cell subsets (ImmuNexUT)", "EAS", best(S1, gene="LYN", ancestry="EAS")),
    ("LYN", "B-cell subsets (ImmuNexUT)", "EUR", best(S1, gene="LYN", ancestry="EUR")),
    ("FCRL5", "B cells (ImmuNexUT)", "EAS", best(S1, gene="FCRL5", ancestry="EAS")),
    ("IRF8", "pDC (ImmuNexUT)", "EUR", best(S1, gene="IRF8", ancestry="EUR")),
    ("TNFSF4", "Naive B (ImmuNexUT)", "EUR", best(S1, gene="TNFSF4", ancestry="EUR")),
    ("LIF", "Monocytes (eQTL Catalogue)", "EAS", best(S4, gene="LIF", ancestry="EAS")),
]
for gene, layer, anc, r in main:
    d = {"Gene or protein": gene, "Layer": layer, "GWAS": anc}
    for p in P12:
        d[f"p12={p:.0e}"] = round(float(rescale(r["PP.H4"], p / 1e-5)), 3)
    rows.append(d)
prior_tab = pd.DataFrame(rows)
out["T1_prior_main"] = prior_tab

# ---------- 2. window counts by prior ----------
windows23 = sorted(set(strata.window))
def count_windows(df, p, thr=0.8):
    pp = rescale(df["PP.H4"], p / 1e-5)
    return df.loc[pp >= thr, "window"].nunique()

rows = []
for p in P12:
    rows.append({
        "p12": f"{p:.0e}",
        "Immune cells (ImmuNexUT), windows of 23": count_windows(S1, p) if p <= 1e-5 else f">={count_windows(S1, p)}",
        "Immune cells excluding CFH": count_windows(S1[S1.gene != "CFH"], p) if p <= 1e-5 else f">={count_windows(S1[S1.gene != 'CFH'], p)}",
        "Kidney (NephQTL2), windows of 24": count_windows(S2, p),
        "Kidney excluding CFH window, of 23": count_windows(S2[S2.window != "CFH"], p),
        "Kidney tests reaching 0.8 (of 4,230)": int((rescale(S2["PP.H4"], p / 1e-5) >= 0.8).sum()),
    })
out["T2_counts_by_prior"] = pd.DataFrame(rows)

# ---------- 3. matched testing, exact enumeration ----------
LINEAGE = {
    "Monocytes": ["CL_Mono", "Int_Mono", "NC_Mono", "CD16p_Mono"],
    "Neutrophils": ["Neu", "LDG"], "Myeloid DC": ["mDC"], "Plasmacytoid DC": ["pDC"],
    "B cells": ["Naive_B", "USM_B", "SM_B", "DN_B"], "Plasmablasts": ["Plasmablast"],
    "CD4 T": ["Naive_CD4", "Mem_CD4", "Th1", "Th2", "Th17", "Tfh"],
    "Treg": ["Fr_I_nTreg", "Fr_II_eTreg", "Fr_III_T"],
    "CD8 T": ["Naive_CD8", "CM_CD8", "EM_CD8", "TEMRA_CD8", "Mem_CD8"], "NK": ["NK"],
}
cells = sorted(S1.cell.unique())
assert sorted(sum(LINEAGE.values(), [])) == cells


def summarize(values, weights):
    v = np.asarray(values, float); w = np.asarray(weights, float); w = w / w.sum()
    o = np.argsort(v); v, w = v[o], w[o]; cw = np.cumsum(w)
    q = lambda a: v[np.searchsorted(cw, a - 1e-12)]
    return {"median": q(0.5), "q25": q(0.25), "q75": q(0.75), "min": v.min(), "max": v.max(),
            "P(count<=1)": float(w[v <= 1].sum())}


def matched(df_pos, label):
    pos = df_pos[df_pos["PP.H4"] >= 0.8]
    pos = pos[pos.window.isin(windows23)]
    win_by_cell = {c: set(pos.loc[pos.cell == c, "window"]) for c in cells}
    res = {}
    # cell-type matched: all C(28,2) pairs equally likely
    vals = [len(win_by_cell[a] | win_by_cell[b]) for a, b in combinations(cells, 2)]
    res["Cell-type matched (exact, 378 pairs)"] = summarize(vals, np.ones(len(vals)))
    # lineage matched: 45 lineage pairs, one cell type uniformly within each
    vals, wts = [], []
    for la, lb in combinations(LINEAGE, 2):
        for a in LINEAGE[la]:
            for b in LINEAGE[lb]:
                vals.append(len(win_by_cell[a] | win_by_cell[b]))
                wts.append(1 / (len(LINEAGE[la]) * len(LINEAGE[lb])))
    res["Lineage matched (exact)"] = summarize(vals, wts)
    # test-count matched: hypergeometric per stratum, Poisson-binomial over windows
    npos = pos.groupby(["window", "ancestry"]).size()
    pwin = []
    for w in windows23:
        p_none = 1.0
        for anc in ("EAS", "EUR"):
            s = strata[(strata.window == w) & (strata.ancestry == anc)].iloc[0]
            N, k = int(s.n_immune), int(min(s.n_kidney, s.n_immune))
            K = int(npos.get((w, anc), 0))
            p_none *= comb(N - K, k) / comb(N, k)
        pwin.append(1 - p_none)
    dist = np.array([1.0])
    for p in pwin:
        dist = np.convolve(dist, [1 - p, p])
    support = np.arange(len(dist)); nz = dist > 1e-15
    res["Test-count matched (exact)"] = summarize(support[nz], dist[nz])
    res["Observed, all 28 cell types"] = {"median": len(set(pos.window))}
    t = pd.DataFrame(res).T
    t.insert(0, "Analysis set", label)
    return t

m_all = matched(S1, "All genes (published analysis)")
m_nocfh = matched(S1[S1.gene != "CFH"], "Excluding immune-cell CFH gene tests")
_w = windows23
windows23 = [w for w in _w if w != "CFH"]
m_noreg = matched(S1[S1.window != "CFH"], "Excluding the whole CFH window (22 windows)")
windows23 = _w
out["T3_matched_exact"] = pd.concat([m_all, m_nocfh, m_noreg])

g = X["S8_genes"]
out["T3b_gene_level"] = pd.DataFrame([
    {"Set": "All", "Combinations": len(g), "Immune PP.H4>=0.8": int((g["immune_best_PP.H4"] >= .8).sum()),
     "Kidney PP.H4>=0.8": int((g["kidney_best_PP.H4"] >= .8).sum())},
    {"Set": "Excluding CFH", "Combinations": int((g.gene != "CFH").sum()),
     "Immune PP.H4>=0.8": int(((g["immune_best_PP.H4"] >= .8) & (g.gene != "CFH")).sum()),
     "Kidney PP.H4>=0.8": int(((g["kidney_best_PP.H4"] >= .8) & (g.gene != "CFH")).sum())}])

# ---------- 4. bootstrap CIs for the ancestry decomposition ----------
rng = np.random.default_rng(20261111)
B = 20000
N_EUR, N_EAS = 503, 504
d = S6.copy()
h = lambda p: 2 * p * (1 - p)
fe = rng.binomial(2 * N_EUR, d.f_eur.values, size=(B, len(d))) / (2 * N_EUR)
fa = rng.binomial(2 * N_EAS, d.f_eas.values, size=(B, len(d))) / (2 * N_EAS)
be = rng.normal(d.b_eur.values, d.se_eur.values, size=(B, len(d)))
ba = rng.normal(d.b_eas.values, d.se_eas.values, size=(B, len(d)))
Ve, Va = h(fe) * be**2, h(fa) * ba**2
dfreq = (h(fa) - h(fe)) * (ba**2 + be**2) / 2
deff = (ba**2 - be**2) * (h(fa) + h(fe)) / 2
ci = lambda a: (np.percentile(a, 2.5), np.percentile(a, 97.5))
nh = ~d.hla.values
rows = []
for i, r in d.iterrows():
    rows.append({"Signal": r.locus, "HLA": bool(r.hla),
                 "V_EUR": r.V_eur, "V_EUR 95% CI": "%.4f-%.4f" % ci(Ve[:, i]),
                 "V_EAS": r.V_eas, "V_EAS 95% CI": "%.4f-%.4f" % ci(Va[:, i]),
                 "Difference": r.dV, "Difference 95% CI": "%.4f to %.4f" % ci(Va[:, i] - Ve[:, i]),
                 "Frequency component 95% CI": "%.4f to %.4f" % ci(dfreq[:, i]),
                 "Effect component 95% CI": "%.4f to %.4f" % ci(deff[:, i])})
out["T4_decomp_per_signal"] = pd.DataFrame(rows)
sE, sA = Ve[:, nh].sum(1), Va[:, nh].sum(1)
sF, sD = dfreq[:, nh].sum(1), deff[:, nh].sum(1)
tot = {
    "Sum V_EUR (25 non-HLA)": (d.V_eur[nh].sum(), *ci(sE)),
    "Sum V_EAS (25 non-HLA)": (d.V_eas[nh].sum(), *ci(sA)),
    "Relative excess in EAS (%)": (100 * (d.V_eas[nh].sum() / d.V_eur[nh].sum() - 1), *ci(100 * (sA / sE - 1))),
    "Frequency share of difference (%)": (100 * d.d_freq[nh].sum() / d.dV[nh].sum(), *ci(100 * sF / (sA - sE))),
    "P(EAS sum > EUR sum)": (np.nan, float((sA > sE).mean()), np.nan),
}
out["T4b_decomp_totals"] = pd.DataFrame(tot, index=["Estimate", "2.5%", "97.5%"]).T

# ---------- 5. Steiger directionality ----------
NGWAS = {"EUR": 5556 + 21178, "EAS": 4590 + 7573}
NQTL = {"eQTL": 1019, "pQTL": 1384}
rows = []
for gene, layer, anc in [("TNFSF13", "eQTL", "EAS"), ("TNFSF13", "pQTL", "EAS"), ("TNFSF13", "eQTL", "EUR"),
                         ("TNFSF13", "pQTL", "EUR"), ("CFH", "pQTL", "EUR"), ("CFH", "pQTL", "EAS"),
                         ("CFHR5", "pQTL", "EUR")]:
    r = best(S5, gene=gene, layer=layer, ancestry=anc)
    z_exp = stats.norm.isf(r.min_qtl_p / 2)
    z_out = r.mr_beta / r.mr_se
    r2e = z_exp**2 / (z_exp**2 + NQTL[layer] - 2)
    r2o = z_out**2 / (z_out**2 + NGWAS[anc] - 2)
    rows.append({"Exposure": f"{gene} {layer}", "GWAS": anc, "Instrument": r.qtl_lead,
                 "r2 exposure": round(r2e, 4), "r2 IgAN": round(r2o, 5),
                 "Direction exposure->IgAN": bool(r2e > r2o),
                 "Wald OR": r.WaldRatio_OR, "PP.H4": r["PP.H4"]})
out["T5_steiger"] = pd.DataFrame(rows)

with pd.ExcelWriter(OUT) as w:
    for k, v in out.items():
        v.to_excel(w, sheet_name=k[:31], index=k in ("T3_matched_exact", "T4b_decomp_totals"))
for k, v in out.items():
    print("=====", k); print(v.to_string())

# component sums (added)
print("diff", d.dV[nh].sum(), ci(sA - sE))
print("freq comp", d.d_freq[nh].sum(), ci(sF))
print("effect comp", d.d_effect[nh].sum(), ci(sD))

# numeric CI table for the Figure 2 script (added)
ci_tab = pd.DataFrame({"locus": d.locus,
    "V_eur_lo": np.percentile(Ve, 2.5, 0), "V_eur_hi": np.percentile(Ve, 97.5, 0),
    "V_eas_lo": np.percentile(Va, 2.5, 0), "V_eas_hi": np.percentile(Va, 97.5, 0),
    "dV_lo": np.percentile(Va - Ve, 2.5, 0), "dV_hi": np.percentile(Va - Ve, 97.5, 0)})
ci_tab.to_csv("Aim1_decomposition_bootstrapCI.tsv", sep="\t", index=False)
