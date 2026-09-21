## ---- Fig. 6: complement risk acts through plasma factor H (needs proj_dir, atl_dir; uses results/figures/Fig4_source_data.tsv from the earlier CFH deletion figure) ----
library(data.table); library(ggplot2); library(patchwork)
FAM <- if (capabilities("aqua")) "Arial" else "sans"   # Kidney International: Arial; lowercase panel labels
fig_dir <- file.path(proj_dir, "results", "figures")
COL <- c(EUR = "#2a78d6", EAS = "#eb6834")
th <- theme_classic(base_size = 7, base_family = FAM) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 6.5, hjust = 0),
        plot.tag = element_text(face = "bold", size = 10), legend.key.size = unit(2.5, "mm"))
w <- "CFH"

## A. CFHR3-CFHR1 deletion: frequency and LD with rs6677604 / kidney CFHR1 eQTL lead, 1000 Genomes super-populations
s4 <- fread(file.path(fig_dir, "Fig4_source_data.tsv"))
s4[, pop := factor(pop, c("EUR", "EAS", "SAS", "AMR", "AFR"))]
s4[, hl := fifelse(pop %in% c("EUR", "EAS"), as.character(pop), "Other")]
pA1 <- ggplot(s4[what == "Deletion frequency"], aes(pop, value, fill = hl)) + geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.4, size = 2) +
  scale_fill_manual(values = c(COL, Other = "#bfbfbf"), guide = "none") +
  scale_y_continuous(limits = c(0, 0.52), expand = c(0, 0)) +
  labs(x = NULL, y = "CFHR3–CFHR1 deletion frequency", tag = "a") + th
ldd <- s4[what != "Deletion frequency"][, marker := fifelse(grepl("GWAS", what), "GWAS lead rs6677604", "Kidney CFHR1 eQTL lead")]
pA2 <- ggplot(ldd, aes(pop, value, shape = marker)) + geom_point(size = 1.8) +
  scale_shape_manual(values = c(16, 2), name = NULL) + scale_y_continuous(limits = c(0, 1.05)) +
  labs(x = NULL, y = expression(r^2 ~ "with deletion")) + th + theme(legend.position = "bottom", legend.direction = "vertical")

## B. PP.H4 at the CFH locus across layers
imx <- fread(file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"))[window == w]
kid <- fread(file.path(atl_dir, "Atlas_kidney_coloc_all.tsv.gz"))[window == w]
ecr <- fread(file.path(atl_dir, "Atlas_eqtlcatalogue_coloc.tsv.gz"))[window == w]
job <- fread(file.path(atl_dir, "Atlas_JOB_pqtl_eqtl.tsv"))[window == w]
best <- function(x, lab) x[, .(PP.H4 = max(PP.H4)), by = ancestry][, lab := lab]
PB <- rbindlist(list(
  best(kid[grepl("^CFH", gene)], "Kidney: CFH / CFHR genes (NephQTL2)"),
  best(ecr[tissue == "kidney cortex" & grepl("^CFH", gene)], "Kidney cortex: CFH / CFHR (GTEx)"),
  best(ecr[tissue == "liver" & gene == "CFH"], "Liver: CFH"),
  best(ecr[tissue == "liver" & gene == "CFHR3"], "Liver: CFHR3"),
  best(ecr[tissue == "liver" & gene == "CFHR4"], "Liver: CFHR4"),
  best(job[layer == "pQTL" & gene == "CFH"], "Plasma factor H"),
  best(job[layer == "pQTL" & gene == "CFHR4"], "Plasma FHR-4"),
  best(job[layer == "pQTL" & gene == "CFHR5"], "Plasma FHR-5"),
  best(imx[gene == "CFH"], "Immune cells: CFH (not interpreted)")), use.names = TRUE)
PB[, grp := fcase(grepl("^Kidney", lab), "Kidney", grepl("^Liver", lab), "Liver", grepl("^Plasma", lab), "Plasma", default = "Immune cells")]
PB[, grp := factor(grp, c("Kidney", "Liver", "Plasma", "Immune cells"))]
PB[, lab := factor(lab, rev(unique(lab)))]
pB <- ggplot(PB, aes(PP.H4, lab, colour = ancestry, shape = ancestry)) +
  annotate("rect", xmin = 0.8, xmax = 1.02, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL, labels = c(EAS = "East Asian GWAS", EUR = "European GWAS"), name = NULL) +
  scale_shape_manual(values = c(EAS = 16, EUR = 17), labels = c(EAS = "East Asian GWAS", EUR = "European GWAS"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 0.8, 1), expand = c(0, 0)) +
  labs(x = "PP.H4", y = NULL, tag = "b") + th +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 6), legend.position = "bottom")

## C. direction-of-effect estimates (Wald ratio at the lead pQTL variant): factor H (primary) and FHR-5 (suggestive)
MR <- job[layer == "pQTL" & gene %in% c("CFH", "CFHR5"),
          .(trait = fifelse(gene == "CFH", "Plasma factor H", "Plasma FHR-5"), ancestry,
            or = exp(mr_beta), lo = exp(mr_beta - 1.96 * mr_se), hi = exp(mr_beta + 1.96 * mr_se), p = mr_p)]
MR[, lab := sprintf("%.2f (%.2f-%.2f)", or, lo, hi)]   # single instrument: P equals the GWAS P, so not shown
pC <- ggplot(MR, aes(or, trait, colour = ancestry)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, orientation = "y", position = position_dodge(width = 0.5), linewidth = 0.4) +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  geom_text(aes(x = max(hi) * 1.1, label = lab), hjust = 0, size = 1.9, position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_colour_manual(values = COL, guide = "none") + scale_x_log10() + coord_cartesian(clip = "off") +
  labs(x = "IgAN odds ratio per unit increase\n(direction-of-effect estimate, Wald ratio; log scale)", y = NULL, tag = "c") + th +
  theme(plot.margin = margin(4, 75, 4, 4))

f6 <- ((pA1 / pA2) | (pB / pC + plot_layout(heights = c(3, 1)))) + plot_layout(widths = c(1, 2))
quartz_ok <- capabilities("aqua")
if (quartz_ok) { quartz(type = "pdf", file = file.path(fig_dir, "Fig6_complement.pdf"), width = 180/25.4, height = 130/25.4); print(f6); dev.off()
} else ggsave(file.path(fig_dir, "Fig6_complement.pdf"), f6, width = 180, height = 130, units = "mm", bg = "white", device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig6_complement.tiff"), f6, width = 180, height = 130, units = "mm", dpi = 600, compression = "lzw", bg = "white")
fwrite(PB, file.path(fig_dir, "Fig6B_source_data.tsv"), sep = "\t"); fwrite(MR, file.path(fig_dir, "Fig6C_source_data.tsv"), sep = "\t")
cat("Fig 6 saved\n"); print(PB); print(MR)
