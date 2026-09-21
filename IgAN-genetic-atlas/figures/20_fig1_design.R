## ---- Fig. 1: study design (atlas version). Standalone: needs only proj_dir. ----
library(ggplot2)
fig_dir <- file.path(proj_dir, "results", "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
box <- function(x, y, w, h, lab, fill = "white", bold = FALSE, size = 2.15)
  list(annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill, colour = "grey20", linewidth = 0.3),
       annotate("text", x = x, y = y, label = lab, size = size, lineheight = 0.95, fontface = if (bold) "bold" else "plain"))
arr <- function(x1, y1, x2, y2, lt = 1) annotate("segment", x = x1, y = y1, xend = x2, yend = y2, linewidth = 0.3, linetype = lt,
                                                 arrow = arrow(length = unit(1.2, "mm"), type = "closed"))
f1 <- ggplot() + coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) + theme_void() +
  box(25, 95, 44, 8, "IgAN GWAS, European ancestry\n5,556 cases / 21,178 controls", fill = "#e8f1fc") +
  box(75, 95, 44, 8, "IgAN GWAS, East Asian ancestry\n4,590 cases / 7,573 controls", fill = "#fdeee7") +
  box(50, 83, 92, 7, "30 independent genome-wide significant signals (Kiryluk et al. 2023)\n24 non-HLA windows (lead variant ±500 kb, merged) + HLA; alleles harmonised, GRCh37") +
  arr(25, 91, 25, 86.5) + arr(75, 91, 75, 86.5) +
  ## three branches
  box(15, 62, 27, 25, "A. Ancestry decomposition\n\nvariance contribution per signal\n(log-odds scale)\nV = 2p(1 − p)β²\n\nEAS − EUR difference split into\nallele-frequency and\neffect-size components\n\nz test of β, FDR across 30 signals\n1000 Genomes; CFHR3–CFHR1 deletion", fill = "#f7f7f5") +
  box(50, 62, 38, 25, "B. Colocalization with four QTL resources\n\n1. Japanese immune cells: ImmuNexUT, 28 cell types, n = 416\n2. European cells and tissues: eQTL Catalogue, 25 datasets\n    (OneK1K, BLUEPRINT, Quach, Schmiedel, GTEx)\n    expression + splicing (LeafCutter); blood, liver, kidney\n3. Japanese whole blood (n = 1,019) and\n    plasma proteins, Olink (n = 1,384) (JCTF)\n4. Kidney: NephQTL2 glomerulus (240), tubulointerstitium (311)", fill = "#f7f7f5") +
  box(85, 62, 27, 25, "C. Sequence-based prediction\n\nAlphaGenome\n24 non-HLA lead variants\n(GRCh38; risk vs other allele)\n\nexpression, chromatin,\nsplice sites and junctions\nacross tissue and\ncell-type tracks", fill = "#f7f7f5") +
  arr(15, 79.5, 15, 74.5) + arr(50, 79.5, 50, 74.5) + arr(85, 79.5, 85, 74.5) +
  box(50, 38, 38, 12, "coloc.abf per gene / splicing cluster, per ancestry and layer\n(≥ 50 shared variants; PP.H4 ≥ 0.8 shared, PP.H3 ≥ 0.8 distinct)\nSuSiE-coloc with ancestry-matched LD (up to 5 signals)\ncis-Mendelian randomization (Wald ratio) for blood and plasma traits", fill = "white") +
  arr(50, 49.5, 50, 44) +
  arr(85, 49.5, 85, 38) + annotate("segment", x = 85, y = 38, xend = 69.2, yend = 38, linewidth = 0.3, linetype = 2,
                                   arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  annotate("text", x = 86, y = 43.5, label = "predictions tested\nagainst splicing QTL", size = 1.9, hjust = 0, fontface = "italic") +
  arr(15, 49.5, 15, 13) + arr(50, 32, 50, 13) +
  box(50, 9, 92, 7, "Cell types, tissues and molecules sharing a genetic signal with each IgAN risk locus, by ancestry", fill = "#f0efec", bold = TRUE, size = 2.4)
ggsave(file.path(fig_dir, "Fig1_design.pdf"), f1, width = 180, height = 130, units = "mm", bg = "white", device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig1_design.tiff"), f1, width = 180, height = 130, units = "mm", dpi = 600, compression = "lzw", bg = "white")
cat("Fig 1 saved\n")
