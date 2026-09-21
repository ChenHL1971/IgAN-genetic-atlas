## ---- Fig. 2: ancestry decomposition of the 30 IgAN signals (needs proj_dir, atl_dir) ----
library(data.table); library(ggplot2); library(patchwork)
fig_dir <- file.path(proj_dir, "results", "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
COL <- c(EUR = "#2a78d6", EAS = "#eb6834"); CMP <- c(Frequency = "#8c8c8c", Effect = "#7a1fa2")
th <- theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.tag = element_text(face = "bold", size = 10), legend.key.size = unit(2.5, "mm"))
d <- fread(file.path(atl_dir, "Aim1_ancestry_decomposition.tsv"))
d[, sig := q_fdr < 0.05]
d[, lab := paste0(locus, fifelse(sig, " *", ""))]
d[, lab := factor(lab, lab[order(dV)])]
tot <- d[hla == FALSE, .(EUR = sum(V_eur), EAS = sum(V_eas), freq = sum(d_freq), eff = sum(d_effect))]
cat(sprintf("non-HLA total: EUR %.3f, EAS %.3f (+%.0f%%); frequency %.0f%%, effect %.0f%%\n", tot$EUR, tot$EAS,
            100 * (tot$EAS / tot$EUR - 1), 100 * tot$freq / (tot$freq + tot$eff), 100 * tot$eff / (tot$freq + tot$eff)))

## A. contribution per signal in each ancestry (dumbbell)
A <- melt(d[, .(lab, hla, EUR = V_eur, EAS = V_eas)], id.vars = c("lab", "hla"), variable.name = "anc", value.name = "V")
pA <- ggplot(A, aes(V, lab)) +
  geom_line(aes(group = lab), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = anc), size = 1.6) +
  scale_colour_manual(values = COL, labels = c(EUR = "European", EAS = "East Asian"), name = NULL) +
  labs(x = expression(atop("Variance contribution (log-odds scale)", 2 * p(1 - p) * beta^2)), y = NULL, tag = "A") + th +
  theme(legend.position = "top", axis.text.y = element_text(face = "italic", size = 6)) + guides(colour = guide_legend(nrow = 2))

## B. EAS - EUR difference split into frequency and effect components
B <- melt(d[, .(lab, Frequency = d_freq, Effect = d_effect)], id.vars = "lab", variable.name = "comp", value.name = "val")
pB <- ggplot(B, aes(val, lab, fill = comp)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(width = 0.7, position = position_stack()) +
  scale_fill_manual(values = CMP, name = NULL, labels = c(Frequency = "Allele-frequency component", Effect = "Effect-size component")) +
  labs(x = "Difference\n(East Asian − European)", y = NULL, tag = "B") + th +
  theme(legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  guides(fill = guide_legend(nrow = 2))

## C. per-allele odds ratios, East Asian vs European
pC <- ggplot(d, aes(or_eur, or_eas)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(shape = hla, colour = sig), size = 1.6) +
  ggrepel::geom_text_repel(data = d[sig | locus %in% c("CARD9", "CD28", "CFH", "LYN")], aes(label = locus),
                           size = 2, fontface = "italic", min.segment.length = 0, segment.size = 0.2, max.overlaps = 20) +
  scale_colour_manual(values = c(`TRUE` = "#d62d20", `FALSE` = "grey40"), labels = c(`TRUE` = "Effect differs (FDR < 0.05)", `FALSE` = "No significant difference"), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "Non-HLA", `TRUE` = "HLA"), name = NULL) +
  coord_equal() + labs(x = "Odds ratio, European GWAS", y = "Odds ratio, East Asian GWAS", tag = "C") + th +
  theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal", legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(), legend.text = element_text(size = 5.5))

f2 <- (pA | pB | pC) + plot_layout(widths = c(1, 0.95, 1.3))
quartz_ok <- capabilities("aqua")
if (quartz_ok) { quartz(type = "pdf", file = file.path(fig_dir, "Fig2_ancestry.pdf"), width = 180/25.4, height = 125/25.4); print(f2); dev.off()
} else ggsave(file.path(fig_dir, "Fig2_ancestry.pdf"), f2, width = 180, height = 125, units = "mm", bg = "white", device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig2_ancestry.tiff"), f2, width = 180, height = 125, units = "mm", dpi = 600, compression = "lzw", bg = "white")
cat("Fig 2 saved\n")
