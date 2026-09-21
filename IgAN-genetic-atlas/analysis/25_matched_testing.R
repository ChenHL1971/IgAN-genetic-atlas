## ---- Step 25: immune-vs-kidney matched-testing sensitivity analysis (Supplemental Fig. S1, Table S8) ----
if (!exists("proj_dir")) proj_dir <- "~/Desktop/IgAN_Genetics_Paper2"
if (!exists("atl_dir")) atl_dir <- file.path(proj_dir, "results", "atlas")
## Needs proj_dir and atl_dir (as in the other scripts). Reads the two colocalization tables written by
## 11b_atlas_immune_coloc.R and 13_atlas_kidney_coloc.R; no other inputs.
for (p in c("data.table", "ggplot2", "patchwork", "openxlsx")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(ggplot2); library(patchwork)
FAM <- if (capabilities("aqua")) "Arial" else "sans"   # Kidney International: Arial; lowercase panel labels
fig_dir <- file.path(proj_dir, "results", "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
N_ITER <- 2000; THR <- 0.8; STRONG <- 1e-5
set.seed(20260921)

imx <- fread(file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"))   # window, ancestry, cell, gene, min_eqtl_p, PP.H4 ...
kid <- fread(file.path(atl_dir, "Atlas_kidney_coloc_all.tsv.gz"))   # window, ancestry, tissue (Glom/Tube), gene, ...

## Directly comparable set: windows tested in both resources (the IGH window has no ImmuNexUT genes)
W <- intersect(unique(imx$window), unique(kid$window))
imx <- imx[window %in% W]; kid <- kid[window %in% W]
obs_imm <- imx[PP.H4 >= THR, uniqueN(window)]
obs_kid <- kid[PP.H4 >= THR, uniqueN(window)]      # glomerulus and tubulointerstitium pooled as one kidney resource
cat(sprintf("Windows tested in both: %d | observed windows with PP.H4 >= %.1f: immune %d, kidney %d\n", length(W), THR, obs_imm, obs_kid))

## (A) Cell-type matched: in each iteration, 2 of the 28 ImmuNexUT cell types drawn with equal probability,
##     without replacement (matching the two kidney compartments); the same 2 cell types are used for all windows.
cells <- sort(unique(imx$cell))
A <- vapply(seq_len(N_ITER), function(i) imx[cell %in% sample(cells, 2) & PP.H4 >= THR, uniqueN(window)], integer(1))

## (B) Test-count matched: within each window x GWAS ancestry, draw (without replacement) as many immune
##     gene-cell-type tests as there are kidney gene-compartment tests (all available if fewer).
kn <- kid[, .(nk = .N), by = .(window, ancestry)]
imx2 <- merge(imx, kn, by = c("window", "ancestry"))
imx2[, grp := .GRP, by = .(window, ancestry)]
idx <- split(seq_len(nrow(imx2)), imx2$grp); nk <- imx2[, nk[1], by = grp][order(grp)]$V1
B <- vapply(seq_len(N_ITER), function(i) {
  s <- unlist(Map(function(ix, n) if (length(ix) <= n) ix else sample(ix, n), idx, nk), use.names = FALSE)
  imx2[s][PP.H4 >= THR, uniqueN(window)] }, integer(1))

## (A2) Lineage-level sensitivity: the 28 ImmuNexUT cell types are not independent (many are subsets of one lineage).
##      Collapse them into 10 broad lineages; in each draw pick 2 of 10 lineages with equal probability (without
##      replacement) and one cell type at random within each, so closely related subsets do not add weight.
##      Also enumerate all 45 lineage pairs using every cell type in both lineages (most generous immune side).
##      Drawn after A and B so the seeded results of A and B are unchanged.
lin_map <- c(CL_Mono = "Monocytes", CD16p_Mono = "Monocytes", Int_Mono = "Monocytes", NC_Mono = "Monocytes",
             Neu = "Neutrophils", LDG = "Neutrophils", mDC = "Myeloid DC", pDC = "Plasmacytoid DC",
             Naive_B = "B cells", USM_B = "B cells", SM_B = "B cells", DN_B = "B cells", Plasmablast = "Plasmablasts",
             Naive_CD4 = "CD4 T", Mem_CD4 = "CD4 T", Th1 = "CD4 T", Th2 = "CD4 T", Th17 = "CD4 T", Tfh = "CD4 T",
             Fr_I_nTreg = "Treg", Fr_II_eTreg = "Treg", Fr_III_T = "Treg",
             Naive_CD8 = "CD8 T", CM_CD8 = "CD8 T", EM_CD8 = "CD8 T", TEMRA_CD8 = "CD8 T", Mem_CD8 = "CD8 T", NK = "NK")
stopifnot(all(cells %in% names(lin_map)))
imx[, lineage := lin_map[cell]]
lins <- sort(unique(lin_map)); by_lin <- split(names(lin_map), lin_map)
A2 <- vapply(seq_len(N_ITER), function(i) {
  ct <- vapply(sample(lins, 2), function(l) { x <- by_lin[[l]]; x[sample.int(length(x), 1)] }, character(1))
  imx[cell %in% ct & PP.H4 >= THR, uniqueN(window)] }, integer(1))
pairs <- combn(lins, 2)
A3 <- apply(pairs, 2, function(p) imx[lineage %in% p & PP.H4 >= THR, uniqueN(window)])
lin_pairs <- data.table(lineage_1 = pairs[1, ], lineage_2 = pairs[2, ], windows_colocalized = A3)[order(windows_colocalized)]
print(lin_pairs[1:5])

## A window counts as colocalized if any sampled test in either GWAS ancestry reaches PP.H4 >= 0.8.
## Empirical one-sided P: probability that a matched immune draw gives a count as low as or lower than kidney,
## with the +1 correction, (k + 1) / (N + 1).
emp <- function(x) (sum(x <= obs_kid) + 1) / (N_ITER + 1)
summ <- rbind(
  data.table(analysis = "Cell-type matched (2 of 28 immune cell types)", median = median(A), q25 = quantile(A, .25), q75 = quantile(A, .75),
             min = min(A), max = max(A), n_le_2 = sum(A <= 2), kidney_observed = obs_kid, empirical_P = emp(A)),
  data.table(analysis = "Test-count matched (per window x ancestry)", median = median(B), q25 = quantile(B, .25), q75 = quantile(B, .75),
             min = min(B), max = max(B), n_le_2 = sum(B <= 2), kidney_observed = obs_kid, empirical_P = emp(B)),
  data.table(analysis = "Lineage matched (2 of 10 lineages, 1 cell type each)", median = median(A2), q25 = quantile(A2, .25), q75 = quantile(A2, .75),
             min = min(A2), max = max(A2), n_le_2 = sum(A2 <= 2), kidney_observed = obs_kid, empirical_P = emp(A2)),
  data.table(analysis = "All 45 lineage pairs (all cell types in both lineages; exhaustive)", median = median(A3), q25 = quantile(A3, .25), q75 = quantile(A3, .75),
             min = min(A3), max = max(A3), n_le_2 = sum(A3 <= 2), kidney_observed = obs_kid, empirical_P = NA_real_))
## exact test-count matching check: strata where immune tests < kidney tests would use all immune tests
strata <- merge(imx[, .(n_immune = .N), by = .(window, ancestry)], kn, by = c("window", "ancestry"))
cat(sprintf("Test-count matching: %d window x ancestry strata; immune tests fewer than kidney tests in %d (immune/kidney ratio %.1f-%.1f)\n",
            nrow(strata), strata[n_immune < nk, .N], min(strata$n_immune / strata$nk), max(strata$n_immune / strata$nk)))
print(summ)

## (C) eQTL detectability: windows with at least one gene with eQTL P < 1e-5 in each resource
det <- data.table(resource = c("ImmuNexUT", "NephQTL2"),
                  windows_with_strong_eqtl = c(imx[min_eqtl_p < STRONG, uniqueN(window)], kid[min_eqtl_p < STRONG, uniqueN(window)]))
print(det)

## (D) Gene-level comparison: same gene, window and GWAS ancestry with eQTL P < 1e-5 in both resources;
##     best PP.H4 across immune cell types and across the two kidney compartments. No further direction or overlap filter
##     beyond the >= 50 shared variants required for any colocalization test.
ig <- imx[min_eqtl_p < STRONG, .(immune_best_PP.H4 = max(PP.H4), immune_best_cell = cell[which.max(PP.H4)]), by = .(window, ancestry, gene)]
kg <- kid[min_eqtl_p < STRONG, .(kidney_best_PP.H4 = max(PP.H4), kidney_best_PP.H3 = PP.H3[which.max(PP.H4)],
                                 kidney_best_compartment = tissue[which.max(PP.H4)]), by = .(window, ancestry, gene)]
D <- merge(ig, kg, by = c("window", "ancestry", "gene"))[order(-immune_best_PP.H4)]
cat(sprintf("Gene-level: %d combinations; immune PP.H4 >= 0.8: %d; kidney PP.H4 >= 0.8: %d\n",
            nrow(D), D[immune_best_PP.H4 >= THR, .N], D[kidney_best_PP.H4 >= THR, .N]))
print(D[immune_best_PP.H4 >= THR])

## ---- Supplemental Fig. S1: null distributions ----
th <- theme_classic(base_size = 7, base_family = FAM) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.tag = element_text(face = "bold", size = 10), plot.title = element_text(size = 7, face = "bold"))
panel <- function(x, ttl, tag) {
  d <- data.table(k = x)[, .N, by = k]
  ggplot(d, aes(k, N)) + geom_col(width = 0.8, fill = "#9aa7b5") +
    geom_vline(xintercept = obs_kid, colour = "#c0392b", linewidth = 0.5) +
    annotate("text", x = obs_kid + 0.3, y = max(d$N) * 1.12, label = sprintf("Kidney tissue\n(observed = %d)", obs_kid),
             hjust = 0, size = 2.1, colour = "#c0392b", lineheight = 0.9) +
    geom_vline(xintercept = median(x), linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = obs_imm, linetype = 3, linewidth = 0.3, colour = "grey30") +
    annotate("text", x = obs_imm + 0.3, y = max(d$N) * 1.12, label = sprintf("All 28 immune-cell\ntypes (observed = %d)", obs_imm),
             hjust = 0, size = 2.1, colour = "grey30", lineheight = 0.9) +
    annotate("text", x = median(x) + 0.3, y = max(d$N) * 1.12, label = sprintf("Median %g\n(IQR %g-%g)", median(x), quantile(x, .25), quantile(x, .75)),
             hjust = 0, size = 2.1, lineheight = 0.9) +
    scale_x_continuous(limits = c(-0.5, length(W) + 0.5), breaks = seq(0, length(W), 5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = sprintf("Windows with PP.H4 ≥ 0.8 (of %d tested in both resources)", length(W)), y = sprintf("Draws (of %d)", N_ITER),
         title = ttl, tag = tag) + th }
fS1 <- panel(A, "Two randomly drawn immune-cell types", "a") / panel(B, "Immune tests matched to kidney test number per window", "b") /
  panel(A2, "Two randomly drawn immune lineages (one cell type each)", "c")
if (capabilities("aqua")) { quartz(type = "pdf", file = file.path(fig_dir, "FigS1_matched_testing.pdf"), width = 120/25.4, height = 160/25.4); print(fS1); dev.off()
} else ggsave(file.path(fig_dir, "FigS1_matched_testing.pdf"), fS1, width = 120, height = 160, units = "mm", bg = "white", device = cairo_pdf)
ggsave(file.path(fig_dir, "FigS1_matched_testing.tiff"), fS1, width = 120, height = 160, units = "mm", dpi = 600, compression = "lzw", bg = "white")

## ---- Supplemental Table S8 ----
S8 <- list(summary = summ, detectability = det, gene_level = D,
           iterations = data.table(iteration = seq_len(N_ITER), cell_type_matched = A, test_count_matched = B, lineage_matched = A2))
fwrite(S8$iterations, file.path(atl_dir, "Matched_testing_iterations.tsv"), sep = "\t")
fwrite(S8$summary, file.path(atl_dir, "Matched_testing_summary.tsv"), sep = "\t")
fwrite(S8$gene_level, file.path(atl_dir, "Matched_testing_gene_level.tsv"), sep = "\t")
fwrite(lin_pairs, file.path(atl_dir, "Matched_testing_lineage_pairs.tsv"), sep = "\t")
fwrite(strata, file.path(atl_dir, "Matched_testing_strata.tsv"), sep = "\t")
cat("Step 25 done\n")
