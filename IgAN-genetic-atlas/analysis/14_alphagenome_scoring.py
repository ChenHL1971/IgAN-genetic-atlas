#!/usr/bin/env python3
"""
Paper 2 (atlas) - Step 14: deep-learning variant-effect prediction with AlphaGenome
(Avsec et al., Nature 2026; https://github.com/google-deepmind/alphagenome).

For each IgAN lead variant (non-HLA signals from Step 10, plus any extra variants in
extra_variants.tsv, e.g. SuSiE credible-set members later) this script:
  1. lifts hg19 -> hg38 and reads the hg38 reference base from the UCSC API
     (so REF/ALT are correct; ALT is then oriented to the IgAN risk allele),
  2. scores the variant with AlphaGenome (1 Mb context) for expression (RNA-seq, CAGE),
     splicing, and chromatin (ATAC, DNase, histone) across all cell types/tissues,
  3. writes every score (risk-allele oriented) and a compact summary by cell group:
     monocyte, dendritic cell, neutrophil, B cell, T cell, NK cell, kidney, liver, intestine.

Setup (once):
  git clone https://github.com/google-deepmind/alphagenome.git && pip install ./alphagenome
  pip install pyliftover pandas requests matplotlib
  API key (free, non-commercial): https://deepmind.google.com/science/alphagenome
  export ALPHAGENOME_API_KEY=...        (or paste it when asked)

usage:  python3 14_alphagenome_scoring.py  PROJECT_DIR  [--only rs3803800]  [--test]
"""
import os, sys, json, time, argparse
import pandas as pd, requests

ap = argparse.ArgumentParser()
ap.add_argument("proj", nargs="?", default=os.path.expanduser("~/Desktop/IgAN_Genetics_Paper2"))
ap.add_argument("--only", default=None, help="comma-separated rsIDs to score (e.g. rs3803800)")
ap.add_argument("--test", action="store_true", help="score rs3803800 only and print a short check")
args = ap.parse_args()
atl = os.path.join(args.proj, "results", "atlas")
out = os.path.join(atl, "alphagenome"); os.makedirs(os.path.join(out, "cache"), exist_ok=True)

# ---------------------------------------------------------------- variants ----
a1 = pd.read_csv(os.path.join(atl, "Aim1_ancestry_decomposition.tsv"), sep="\t")
v = a1.loc[~a1["hla"].astype(bool), ["locus", "rsid", "chr", "pos", "risk", "other"]].copy()
xf = os.path.join(atl, "extra_variants.tsv")          # optional: locus rsid chr pos risk other (hg19)
if os.path.exists(xf):
    v = pd.concat([v, pd.read_csv(xf, sep="\t")], ignore_index=True)
if args.test: args.only = "rs3803800"
if args.only: v = v[v["rsid"].isin(args.only.split(","))]
v = v.dropna(subset=["risk", "other"]).drop_duplicates("rsid").reset_index(drop=True)
print(f"{len(v)} variants to score")

from pyliftover import LiftOver
lo = LiftOver("hg19", "hg38")
def ref38(chrom, pos):                                  # 1-based position -> base on hg38
    r = requests.get("https://api.genome.ucsc.edu/getData/sequence",
                     params={"genome": "hg38", "chrom": chrom, "start": pos - 1, "end": pos}, timeout=60)
    r.raise_for_status(); return r.json()["dna"].upper()
rows = []
for _, r in v.iterrows():
    m = lo.convert_coordinate(f"chr{int(r.chr)}", int(r.pos) - 1)
    if not m: print(f"  {r.rsid}: liftOver failed, skipped"); continue
    c38, p38 = m[0][0], m[0][1] + 1
    base = ref38(c38, p38); risk, oth = str(r.risk).upper(), str(r.other).upper()
    if len(risk) != 1 or len(oth) != 1: print(f"  {r.rsid}: indel, skipped"); continue
    if base == oth:   ref, alt, sign = oth, risk, 1       # ALT = risk allele
    elif base == risk: ref, alt, sign = risk, oth, -1     # ALT = protective -> flip sign
    else: print(f"  {r.rsid}: hg38 base {base} matches neither {risk}/{oth} (strand?), skipped"); continue
    rows.append(dict(locus=r.locus, rsid=r.rsid, chrom=c38, pos38=p38, ref=ref, alt=alt, sign=sign, risk=risk))
v38 = pd.DataFrame(rows); v38.to_csv(os.path.join(out, "variants_hg38.tsv"), sep="\t", index=False)
print(v38.to_string(index=False))

# --------------------------------------------------------------- AlphaGenome ----
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
key = os.environ.get("ALPHAGENOME_API_KEY") or input("AlphaGenome API key: ").strip()
model = dna_client.create(key)
L = dna_client.SUPPORTED_SEQUENCE_LENGTHS["SEQUENCE_LENGTH_1MB"]
WANT = ("RNA_SEQ", "CAGE", "SPLICE", "ATAC", "DNASE", "CHIP_HISTONE")
rec = variant_scorers.RECOMMENDED_VARIANT_SCORERS
scorers = [rec[k] for k in rec if any(w in k.upper() for w in WANT)]
print("scorers:", [k for k in rec if any(w in k.upper() for w in WANT)])

frames = []
for _, r in v38.iterrows():
    cf = os.path.join(out, "cache", f"{r.rsid}.tsv.gz")
    if os.path.exists(cf):
        frames.append(pd.read_csv(cf, sep="\t")); print(f"  {r.rsid}: cached"); continue
    var = genome.Variant(chromosome=r.chrom, position=int(r.pos38), reference_bases=r.ref,
                         alternate_bases=r.alt, name=r.rsid)
    t0 = time.time()
    for attempt in range(3):
        try:
            sc = model.score_variant(interval=var.reference_interval.resize(L), variant=var,
                                     variant_scorers=scorers, organism=dna_client.Organism.HOMO_SAPIENS)
            break
        except Exception as e:
            print(f"  {r.rsid}: attempt {attempt + 1} failed ({e}); retrying in 30 s"); time.sleep(30)
    else:
        print(f"  {r.rsid}: gave up"); continue
    d = variant_scorers.tidy_scores([sc])
    d.insert(0, "rsid", r.rsid); d.insert(1, "locus", r.locus)
    d["risk_raw"] = d["raw_score"] * r.sign           # risk-allele oriented
    d["risk_quantile"] = d["quantile_score"] * r.sign
    d.to_csv(cf, sep="\t", index=False); frames.append(d)
    print(f"  {r.rsid}: {len(d)} scores in {time.time() - t0:.0f} s")
if not frames: sys.exit("no scores")
allsc = pd.concat(frames, ignore_index=True)
allsc.to_csv(os.path.join(out, "AlphaGenome_all_scores.tsv.gz"), sep="\t", index=False)

# ------------------------------------------------------------------ summary ----
GROUPS = {  # regex on biosample_name / gtex_tissue (case-insensitive)
    "Monocyte": r"monocyte", "Dendritic cell": r"dendritic", "Neutrophil": r"neutrophil",
    "B cell": r"\bB cell|B-cell|naive B|memory B", "T cell": r"T cell|T-cell|T-helper|regulatory T|CD4|CD8",
    "NK cell": r"natural killer|NK cell", "Kidney": r"kidney|renal|glomerul|podocyte|proximal tubul",
    "Liver": r"liver|hepatocyte", "Intestine": r"colon|intestin|ileum|small bowel|duoden|rectum"}
lab = allsc["biosample_name"].fillna("") + " " + allsc.get("gtex_tissue", pd.Series("", index=allsc.index)).fillna("")
allsc["cell_group"] = None
for g, rx in GROUPS.items():
    allsc.loc[allsc["cell_group"].isna() & lab.str.contains(rx, case=False, regex=True), "cell_group"] = g
allsc["modality"] = allsc["output_type"].astype(str).str.upper()
sub = allsc.dropna(subset=["cell_group"])
def top(d):  # strongest risk-allele effect in the group
    i = d["risk_quantile"].abs().idxmax(); return d.loc[i, ["gene_name", "risk_quantile", "risk_raw", "biosample_name", "track_name"]]
summ = (sub.groupby(["locus", "rsid", "modality", "cell_group"]).apply(top).reset_index())
summ.to_csv(os.path.join(out, "AlphaGenome_summary_by_cell_group.tsv"), sep="\t", index=False)

# expression heatmap (RNA-seq): strongest gene effect per variant x cell group
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
e = summ[summ["modality"].str.contains("RNA")]
if len(e):
    hm = e.pivot_table(index="locus", columns="cell_group", values="risk_quantile", aggfunc="first")
    hm = hm[[g for g in GROUPS if g in hm.columns]]
    fig, ax = plt.subplots(figsize=(1.2 + 0.55 * hm.shape[1], 0.9 + 0.28 * hm.shape[0]))
    im = ax.imshow(hm.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(hm.shape[1])); ax.set_xticklabels(hm.columns, rotation=50, ha="right", fontsize=7)
    ax.set_yticks(range(hm.shape[0])); ax.set_yticklabels(hm.index, fontsize=7, style="italic")
    for i in range(hm.shape[0]):
        for j in range(hm.shape[1]):
            x = hm.values[i, j]
            if pd.notna(x) and abs(x) >= 0.9: ax.text(j, i, f"{x:.2f}", ha="center", va="center", fontsize=5.5, color="white")
    cb = fig.colorbar(im, ax=ax, shrink=0.6); cb.set_label("Risk-allele effect on expression\n(AlphaGenome quantile score)", fontsize=7)
    ax.set_title("Predicted expression effects of IgAN lead variants", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "AlphaGenome_expression_heatmap.pdf")); plt.close(fig)

if args.test:
    t = allsc[(allsc["gene_name"] == "TNFSF13") & allsc["cell_group"].notna()]
    print("\n[test] rs3803800 -> TNFSF13, strongest effects by cell group:")
    print(t.loc[t.groupby(["modality", "cell_group"])["risk_quantile"].apply(lambda s: s.abs().idxmax())]
          [["modality", "cell_group", "biosample_name", "risk_quantile"]].to_string(index=False))
print("Step 14 done:", out)
