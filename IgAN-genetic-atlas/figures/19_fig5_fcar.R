## ---- Fig. 5: FCAR splice-site variant (paste after Step 11b Part 1; needs get_gwas, kg_window, lift19, panel, atl_dir) ----
library(data.table); library(ggplot2); library(patchwork)
fig_dir <- file.path(proj_dir, "results", "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
COL <- c(EUR = "#2a78d6", EAS = "#eb6834")
LDC <- c("#bfbfbf", "#86b6ef", "#3fae5a", "#f0a030", "#d62d20")
th <- theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 6.5, hjust = 0),
        plot.tag = element_text(face = "bold", size = 10), legend.key.size = unit(2.5, "mm"))
w <- "FCAR"; lead19 <- "19:55397217"; pos38 <- 54885805; jx <- "19:54885525:54885805"   # rs1865097; junction (GRCh38)
ec_dir <- file.path(atl_dir, "eqtlcat")

## A. regional: East Asian GWAS and the junction sQTL (Quach monocytes, GTEx blood), LD with rs1865097 in EAS
k <- kg_window(w); eas <- intersect(panel[super_pop == "EAS", sample], colnames(k$D)); X <- k$D[, eas, drop = FALSE]
ld <- data.table(SNP = k$meta$SNP, r2 = suppressWarnings(as.vector(cor(t(X), X[match(lead19, k$meta$SNP), ])^2)))
bin_r2 <- function(r) cut(r, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">=0.8"))
g <- get_gwas(w, "EAS")[, .(SNP, pos = as.integer(sub(".*:", "", SNP)), P, track = "IgAN GWAS, East Asian")]
dd <- data.table(did = c("QTD000413", "QTD000360"), track = c("Junction sQTL, monocytes (Quach)", "Junction sQTL, whole blood (GTEx)"))
sq <- rbindlist(lapply(seq_len(nrow(dd)), function(i) {
  x <- readRDS(file.path(ec_dir, sprintf("%s_%s.rds", dd$did[i], w)))[startsWith(molecular_trait_id, paste0(jx, ":"))]
  x[, .(position, P = pvalue, beta, se, did = dd$did[i], track = dd$track[i])] }))
up <- unique(sq[, .(position)])[, pos := lift19("chr19", position)]
sq <- merge(sq, up[!is.na(pos)], by = "position")[, SNP := paste0("19:", pos)]
R <- rbind(g, sq[, .(SNP, pos, P, track)])[P > 0]
R <- merge(R, ld, by = "SNP", all.x = TRUE)[, r2b := bin_r2(r2)][is.na(r2b), r2b := "<0.2"]
R[, track := factor(track, c("IgAN GWAS, East Asian", "Junction sQTL, monocytes (Quach)", "Junction sQTL, whole blood (GTEx)"))]
R <- R[abs(pos - 55397217) <= 250000][order(r2b)]
pA <- ggplot(R, aes(pos / 1e6, -log10(P))) +
  geom_point(aes(fill = r2b), shape = 21, size = 0.9, stroke = 0.1, colour = "grey30") +
  geom_point(data = R[SNP == lead19], shape = 23, size = 2, fill = "#7a1fa2", colour = "black", stroke = 0.3) +
  facet_wrap(~track, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = setNames(LDC, levels(R$r2b)), name = expression(r^2 ~ "(EAS)"), drop = FALSE) +
  labs(x = "Chromosome 19 position (Mb, GRCh37)", y = expression(-log[10] ~ italic(P)), tag = "A") + th +
  theme(legend.position = "bottom", legend.title = element_text(size = 6)) + guides(fill = guide_legend(nrow = 2, override.aes = list(size = 1.8)))

## B. AlphaGenome: predicted change in use of this junction, top 15 tracks (unsigned)
ag <- fread(file.path(atl_dir, "alphagenome", "AlphaGenome_all_scores.tsv.gz"))
a <- ag[rsid == "rs1865097" & gene_name == "FCAR" & grepl("junction", variant_scorer, ignore.case = TRUE) &
        junction_Start == 54885525][order(-abs(raw_score))]
a <- a[, .SD[1], by = biosample_name][order(-abs(raw_score))][1:min(.N, 15)]
a[, hl := fifelse(grepl("monocyte|blood|spleen", biosample_name, ignore.case = TRUE), "Myeloid / blood", "Other")]
a[, biosample_name := factor(biosample_name, rev(biosample_name))]
pB <- ggplot(a, aes(abs(raw_score), biosample_name, fill = hl)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c(`Myeloid / blood` = "#eb6834", Other = "#bfbfbf"), name = NULL) +
  scale_x_continuous(expand = expansion(c(0, 0.05))) +
  labs(x = "Absolute change in junction use (AlphaGenome)", y = NULL, tag = "B") + th + theme(legend.position = "bottom")

## C. sQTL effect of the A (risk) allele on FCAR junctions at rs1865097
jj <- rbindlist(lapply(c("QTD000413", "QTD000360"), function(d)
  readRDS(file.path(ec_dir, sprintf("%s_%s.rds", d, w)))[position == pos38 & gene_id == "ENSG00000186431"][, did := d]))
jj[, junction := sub("^19:(\\d+):(\\d+):.*$", "\\1-\\2", molecular_trait_id)]
jj[, beta_A := fifelse(alt == "A", beta, -beta)]      # effect per A (risk) allele
jj[, study := fifelse(did == "QTD000413", "Monocytes (Quach)", "Whole blood (GTEx)")]
jj[, target := fifelse(startsWith(molecular_trait_id, paste0(jx, ":")), "Affected junction", "Other FCAR junctions")]
pC <- ggplot(jj, aes(beta_A, junction, colour = target)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = beta_A - 1.96 * se, xmax = beta_A + 1.96 * se), width = 0.2, orientation = "y", linewidth = 0.4) +
  geom_point(size = 1.5) + facet_wrap(~study, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(`Affected junction` = "#d62d20", `Other FCAR junctions` = "grey45"), name = NULL) +
  labs(x = "Splicing QTL effect of A allele (β, 95% CI)", y = "Junction (GRCh38)", tag = "C") + th + theme(legend.position = "bottom")

## D. PP.H4 at the FCAR locus across layers
imx <- fread(file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"))[window == w & gene == "FCAR"]
ecr <- fread(file.path(atl_dir, "Atlas_eqtlcatalogue_coloc.tsv.gz"))[window == w & gene == "FCAR"]
job <- fread(file.path(atl_dir, "Atlas_JOB_pqtl_eqtl.tsv"))[window == w & (gene == "FCAR" | grepl("^KIR", gene))]
kid <- fread(file.path(atl_dir, "Atlas_kidney_coloc_all.tsv.gz"))[window == w & gene == "FCAR"]
PD <- rbindlist(list(
  imx[, .(lab = "FCAR expression, Japanese immune cells", PP.H4 = max(PP.H4)), by = ancestry],
  ecr[quant == "ge", .(lab = "FCAR expression, European cells", PP.H4 = max(PP.H4)), by = ancestry],
  ecr[quant == "leafcutter", .(lab = "FCAR splicing, monocytes/blood", PP.H4 = max(PP.H4)), by = ancestry],
  job[gene == "FCAR" & layer == "eQTL", .(lab = "FCAR expression, whole blood (JCTF)", PP.H4 = max(PP.H4)), by = ancestry],
  job[gene == "FCAR" & layer == "pQTL", .(lab = "Plasma CD89 (FCAR)", PP.H4 = max(PP.H4)), by = ancestry],
  job[grepl("^KIR", gene) & layer == "pQTL", .(lab = "Plasma KIR proteins", PP.H4 = max(PP.H4)), by = ancestry],
  kid[, .(lab = "FCAR expression, kidney", PP.H4 = max(PP.H4)), by = ancestry]), use.names = TRUE)
PD[, lab := factor(lab, rev(unique(lab)))]
pD <- ggplot(PD, aes(PP.H4, lab, colour = ancestry, shape = ancestry)) +
  annotate("rect", xmin = 0.8, xmax = 1.02, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = COL, labels = c(EAS = "East Asian GWAS", EUR = "European GWAS"), name = NULL) +
  scale_shape_manual(values = c(EAS = 16, EUR = 17), labels = c(EAS = "East Asian GWAS", EUR = "European GWAS"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 0.8, 1), expand = c(0, 0)) +
  labs(x = "PP.H4", y = NULL, tag = "D") + th + theme(legend.position = "bottom")

f5 <- wrap_elements(full = (pA | (pB / pC))) / wrap_elements(full = pD) + plot_layout(heights = c(3, 1.1))
ggsave(file.path(fig_dir, "Fig5_FCAR_splice.pdf"), f5, width = 180, height = 170, units = "mm", bg = "white", device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig5_FCAR_splice.tiff"), f5, width = 180, height = 170, units = "mm", dpi = 600, compression = "lzw", bg = "white")
fwrite(jj, file.path(fig_dir, "Fig5C_source_data.tsv"), sep = "\t"); fwrite(PD, file.path(fig_dir, "Fig5D_source_data.tsv"), sep = "\t")
cat("Fig 5 saved\n"); print(jj[, .(study, junction, ref, alt, beta_A = round(beta_A, 3), pvalue = signif(pvalue, 2))]); print(PD)
