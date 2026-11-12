"""Supplementary Figure S2 (usage: python 29_figS2_fcar_transcript.py FigS2_FCAR_transcript.tif): FCAR transcript structure and the rs1865097 junction.
Coordinates: Ensembl canonical transcript ENST00000355524 (FCAR-203), GRCh38."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Arc
import sys

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
OUT = sys.argv[1]

exons = [(54874235, 54874323, "1", "5′UTR / S1"), (54875330, 54875365, "2", "S2"),
         (54885235, 54885525, "3", "EC1"), (54888007, 54888294, "4", "EC2"),
         (54889649, 54891420, "5", "Stalk / TM /\ncytoplasmic / 3′UTR")]
cryptic = (54885805, 54885857)  # cryptic exon start (rs1865097) to in-frame stop codon
start, end = 54873900, 54891700

fig = plt.figure(figsize=(7.2, 6.8))
ax = fig.add_axes([0.04, 0.55, 0.92, 0.42])
ax.set_xlim(start, end); ax.set_ylim(-3.2, 4.4); ax.axis("off")
x = lambda p: p
grey, blue, red, orange = "#555555", "#2F6DB5", "#C0392B", "#E07B39"

ax.text(start, 4.3, "a", fontsize=13, fontweight="bold", va="top")
ax.text(start + 450, 4.25, "FCAR (canonical transcript ENST00000355524; chr19, + strand; GRCh38)", va="top")
ax.plot([exons[0][0], exons[-1][1]], [0, 0], color=grey, lw=1)
for s, e, n, lab in exons:
    w = max(e - s, 160)
    ax.add_patch(Rectangle((s, -0.45), w, 0.9, color=blue, zorder=3))
    ax.text(s + w / 2, -0.9, f"Exon {n}", ha="center", va="top", fontsize=8)
    ax.text(s + w / 2, -1.55, lab, ha="center", va="top", fontsize=7.5, color="#333333")
# cryptic exon
cs, ce = cryptic
ax.add_patch(Rectangle((cs, -0.3), 160, 0.6, color=orange, zorder=3))
ax.annotate("Cryptic\nexon", xy=(cs + 80, -0.3), xytext=(54886750, -0.9), ha="center", va="top", fontsize=7.5, color=orange, arrowprops=dict(arrowstyle="-", color=orange, lw=0.8))

def arc(a, b, h, color, ls="-"):
    mid, wid = (a + b) / 2, (b - a)
    ax.add_patch(Arc((mid, 0.45), wid, h, theta1=0, theta2=180, color=color, lw=1.6, ls=ls))

e3_end = 54885525
arc(e3_end, 54888007, 4.6, blue)
ax.text(54888300, 2.2, "Canonical junction\nexon 3 → exon 4", ha="left", color=blue, fontsize=8)
arc(e3_end, cs, 2.2, red)
ax.annotate("Affected junction\n19:54885525–54885805\n(reduced by risk allele A)", xy=(cs - 120, 1.25),
            xytext=(54879300, 1.75), color=red, fontsize=8, ha="center",
            arrowprops=dict(arrowstyle="-", color=red, lw=0.8))

# zoom panel b: sequence at cryptic acceptor
bx = fig.add_axes([0.04, 0.29, 0.92, 0.24]); bx.axis("off"); bx.set_xlim(0, 100); bx.set_ylim(-0.8, 10)
bx.text(0, 9.8, "b", fontsize=13, fontweight="bold", va="top")
bx.text(4, 9.9, "Intron 3, ~280 bp downstream of exon 3 (reference sequence, GRCh38)", va="top")
seq_intron = "…ttttacctgatagtctttctctgttccag"
fig.canvas.draw(); rend = fig.canvas.get_renderer()
def chain(x0, parts, y):
    xcur = x0
    for txt, col, wt in parts:
        t = bx.text(xcur, y, txt, family="DejaVu Sans Mono", fontsize=8.5, color=col, fontweight=wt, va="center")
        bb = t.get_window_extent(renderer=rend).transformed(bx.transData.inverted())
        if txt == "G": gx = (bb.x0 + bb.x1) / 2
        xcur = bb.x1
    return gx
gx = chain(4, [(seq_intron, grey, "normal"), ("G", red, "bold"),
               ("ATCCCATTCAAGATTTCACATTGCGGCTGGGAGTGGTGGCTCACGCCTGTAA", orange, "normal")], 4.6)
bx.annotate("rs1865097 (19:54885805): G protective / A risk;\nfirst base of the cryptic exon", xy=(gx, 5.3),
            xytext=(gx, 6.6), ha="center", fontsize=7.5, color=red, arrowprops=dict(arrowstyle="-", color=red, lw=0.8))
bx.text(4, 2.0, "AG acceptor at 19:54885803–54885804 after a polypyrimidine tract.\nIn-frame stop codon (TAA) at nucleotides 51–53 of the cryptic exon.", fontsize=7.5, color=grey)
bx.text(4, 0.0, "Reading frame continued from exon 3: …EC1–G-S-H-S-R-F-H-I-A-A-G-S-G-G-S-R-L–stop", fontsize=7.5)

# panel c: predicted products
cx = fig.add_axes([0.04, 0.02, 0.92, 0.26]); cx.axis("off"); cx.set_xlim(0, 100); cx.set_ylim(0, 10)
cx.text(0, 9.8, "c", fontsize=13, fontweight="bold", va="top")
cx.text(4, 9.6, "Predicted protein products", va="top")
def bar(y, segs, label, note):
    xx = 22
    for w, txt, col in segs:
        cx.add_patch(FancyBboxPatch((xx, y - 0.7), w, 1.4, boxstyle="round,pad=0.02,rounding_size=0.3",
                                    color=col, lw=0))
        cx.text(xx + w / 2, y, txt, ha="center", va="center", color="white", fontsize=7.5)
        xx += w + 0.5
    cx.text(4, y, label, va="center", fontsize=8)
    cx.text(xx + 1, y, note, va="center", fontsize=7.5, color="#333333")
bar(6.4, [(5, "SP", "#7F8C8D"), (16, "EC1", blue), (16, "EC2", blue), (6, "Stalk", "#5D8AA8"),
          (5, "TM", "#1B2631"), (8, "Cyto", "#5D8AA8")], "Exon 3 → exon 4\n(canonical)",
    "Membrane CD89\n(287 aa)")
bar(2.8, [(5, "SP", "#7F8C8D"), (16, "EC1", blue), (3, "", orange)], "Exon 3 → cryptic\nexon",
    "137 aa; no EC2, TM or cytoplasmic tail;\nmay be degraded by NMD (not measured here)")
fig.savefig(OUT + ".png", dpi=200)
fig.savefig(OUT + ".tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
