## ---- Fig. 4: the APRIL chain (paste after Step 11b Part 1 + Part 2; needs get_gwas, kg_window, lift19, panel, co_dir, atl_dir) ----
for (p in c("data.table", "ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(ggplot2); library(patchwork)
fig_dir <- file.path(proj_dir, "results", "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
COL <- c(EUR = "#2a78d6", EAS = "#eb6834")
LDC <- c("#bfbfbf", "#86b6ef", "#3fae5a", "#f0a030", "#d62d20")          # r2 bins <0.2 ... >=0.8
w <- "TNFSF12_13"; lead <- "17:7462969"; xr <- c(7.25, 7.70)             # Mb, GRCh37
th <- theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 6.5, hjust = 0),
        plot.tag = element_text(face = "bold", size = 10), legend.key.size = unit(2.5, "mm"))

## LD (r2) with rs3803800 in 1000 Genomes East Asians
k <- kg_window(w); eas <- intersect(panel[super_pop == "EAS", sample], colnames(k$D))
X <- k$D[, eas, drop = FALSE]; li <- match(lead, k$meta$SNP)
ld <- data.table(SNP = k$meta$SNP, r2 = suppressWarnings(as.vector(cor(t(X), X[li, ])^2)))
bin_r2 <- function(r) cut(r, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">=0.8"))

## A. regional tracks: -log10 P vs position
tr <- list()
for (a in c("EAS", "EUR")) { g <- get_gwas(w, a)
  tr[[a]] <- g[, .(SNP, pos = as.integer(sub(".*:", "", SNP)), P, track = sprintf("IgAN GWAS, %s", fifelse(a == "EAS", "East Asian", "European")))] }
e <- fread(file.path(co_dir, paste0("imx_", w, "_hg38.tsv")))[cell == "CL_Mono" & Gene_name == "TNFSF13"]
up <- unique(e[, .(Variant_CHR, Variant_position_start)])[, pos := lift19(Variant_CHR, Variant_position_start)]
e <- merge(e, up[!is.na(pos)], by = c("Variant_CHR", "Variant_position_start"))
tr$imx <- e[, .(SNP = paste(sub("^chr", "", Variant_CHR), pos, sep = ":"), pos, P = nominal_P_value, track = "TNFSF13 eQTL, monocytes")]
jb_out <- file.path(atl_dir, "job")
for (kk in c("eQTL", "pQTL")) {
  j <- fread(file.path(jb_out, paste0(kk, "_windows.tsv.gz")))[gene_name == "TNFSF13"]
  j[, c("c", "p") := tstrsplit(variant_id_hg19, ":", fixed = TRUE, keep = 1:2)]
  tr[[kk]] <- j[, .(SNP = paste(sub("^chr", "", c), p, sep = ":"), pos = as.integer(p), P = pval_nominal,
                    track = fifelse(kk == "eQTL", "TNFSF13 eQTL, whole blood", "Plasma APRIL pQTL"))]
}
R <- rbindlist(tr)[P > 0]
R <- merge(R, ld, by = "SNP", all.x = TRUE)[, r2b := bin_r2(r2)]
R[is.na(r2b), r2b := "<0.2"]
R[, track := factor(track, c("IgAN GWAS, East Asian", "IgAN GWAS, European", "TNFSF13 eQTL, monocytes",
                               "TNFSF13 eQTL, whole blood", "Plasma APRIL pQTL"))]
R <- R[pos / 1e6 >= xr[1] & pos / 1e6 <= xr[2]][order(r2b)]
pA <- ggplot(R, aes(pos / 1e6, -log10(P))) +
  geom_point(aes(fill = r2b), shape = 21, size = 0.9, stroke = 0.1, colour = "grey30") +
  geom_point(data = R[SNP == lead], shape = 23, size = 2, fill = "#7a1fa2", colour = "black", stroke = 0.3) +
  facet_wrap(~track, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_fill_manual(values = setNames(LDC, levels(R$r2b)), name = expression(r^2 ~ "(EAS)"), drop = FALSE) +
  labs(x = "Chromosome 17 position (Mb, GRCh37)", y = expression(-log[10] ~ italic(P)), tag = "A") + th +
  theme(legend.position = "bottom", legend.title = element_text(size = 6)) + guides(fill = guide_legend(nrow = 2, override.aes = list(size = 1.8)))

## B. PP.H4 for TNFSF13 across layers and ancestries
imx <- fread(file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"))[window == w & gene == "TNFSF13"]
lin <- c(CL_Mono = "Monocytes", Int_Mono = "Monocytes", NC_Mono = "Monocytes", CD16p_Mono = "Monocytes", mDC = "Myeloid DC",
         pDC = "Plasmacytoid DC", Neu = "Neutrophils", LDG = "Neutrophils", Naive_B = "B cells", USM_B = "B cells", SM_B = "B cells",
         DN_B = "B cells", Plasmablast = "B cells", NK = "NK cells")
imx[, lab := fifelse(cell %in% names(lin), lin[cell], fifelse(grepl("CD8", cell), "CD8 T cells", "CD4 T cells"))]
b1 <- imx[, .(PP.H4 = max(PP.H4)), by = .(lab, ancestry)][, layer := "Japanese immune cells"]
ecr <- fread(file.path(atl_dir, "Atlas_eqtlcatalogue_coloc.tsv.gz"))[window == w & gene == "TNFSF13"]
b2 <- ecr[, .(PP.H4 = max(PP.H4)), by = .(lab = fifelse(tissue %in% c("monocyte", "CD16+ monocyte"), "Monocytes",
                                                   fifelse(tissue == "blood", "Whole blood", "Other cells")), ancestry)][, layer := "European cells"]
job <- fread(file.path(atl_dir, "Atlas_JOB_pqtl_eqtl.tsv"))[window == w & gene == "TNFSF13"]
b3 <- job[, .(lab = fifelse(layer == "pQTL", "Plasma APRIL", "Whole blood"), ancestry, PP.H4, layer = "Japanese blood")]
kid <- fread(file.path(atl_dir, "Atlas_kidney_coloc_all.tsv.gz"))[window == w & gene == "TNFSF13"]
b4 <- kid[, .(lab = fifelse(tissue == "Glom", "Glomerulus", "Tubulointerstitium"), ancestry, PP.H4, layer = "Kidney")]
PB <- rbindlist(list(b1, b2, b3, b4), use.names = TRUE)
PB[, layer := factor(layer, c("Japanese immune cells", "European cells", "Japanese blood", "Kidney"))]
lab_order <- c("Monocytes", "Myeloid DC", "Plasmacytoid DC", "Neutrophils", "B cells", "CD4 T cells", "CD8 T cells", "NK cells",
               "Other cells", "Whole blood", "Plasma APRIL", "Glomerulus", "Tubulointerstitium")
PB[, lab := factor(lab, rev(lab_order))]
pB <- ggplot(PB, aes(PP.H4, lab, colour = ancestry, shape = ancestry)) +
  annotate("rect", xmin = 0.8, xmax = 1.02, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL, labels = c(EAS = "East Asian GWAS", EUR = "European GWAS"), name = NULL) +
  scale_shape_manual(values = c(EAS = 16, EUR = 17), labels = c(EAS = "East Asian GWAS", EUR = "European GWAS"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 0.8, 1), expand = c(0, 0)) +
  labs(x = "PP.H4 (TNFSF13 / APRIL)", y = NULL, tag = "B") + th +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 6), legend.position = "bottom")

## C. cis-MR at the lead QTL variant
MR <- job[, .(trait = fifelse(layer == "pQTL", "Plasma APRIL", "Whole-blood TNFSF13"), ancestry,
              or = exp(mr_beta), lo = exp(mr_beta - 1.96 * mr_se), hi = exp(mr_beta + 1.96 * mr_se), p = mr_p)]
MR[, lab := sprintf("%.2f (%.2f-%.2f); P = %s", or, lo, hi, formatC(p, format = "e", digits = 1))]
pC <- ggplot(MR, aes(or, trait, colour = ancestry)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, orientation = "y", position = position_dodge(width = 0.5), linewidth = 0.4) +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  geom_text(aes(x = max(hi) * 1.08, label = lab), hjust = 0, size = 1.9, position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_colour_manual(values = COL, guide = "none") + scale_x_log10() + coord_cartesian(clip = "off") +
  labs(x = "IgAN odds ratio per unit increase (cis-MR, log scale)", y = NULL, tag = "C") + th +
  theme(plot.margin = margin(4, 70, 4, 4))

f4 <- pA + (pB / pC + plot_layout(heights = c(3.2, 1))) + plot_layout(widths = c(1.25, 1))
ggsave(file.path(fig_dir, "Fig4_APRIL_chain.pdf"), f4, width = 180, height = 150, units = "mm", bg = "white")
ggsave(file.path(fig_dir, "Fig4_APRIL_chain.tiff"), f4, width = 180, height = 150, units = "mm", dpi = 600, compression = "lzw", bg = "white")
fwrite(PB, file.path(fig_dir, "Fig4B_source_data.tsv"), sep = "\t"); fwrite(MR, file.path(fig_dir, "Fig4C_source_data.tsv"), sep = "\t")
cat("Fig 4 saved\n"); print(R[SNP == lead, .(track, P, r2)]); print(MR)
