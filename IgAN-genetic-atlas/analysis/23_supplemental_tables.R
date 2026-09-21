## ---- Supplemental Tables S1-S8 (one Excel workbook) — needs proj_dir, atl_dir ----
if (!exists("proj_dir")) proj_dir <- "~/Desktop/IgAN_Genetics_Paper2"
if (!exists("atl_dir")) atl_dir <- file.path(proj_dir, "results", "atlas")
for (p in c("data.table", "openxlsx")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(openxlsx)
sub_dir <- file.path(proj_dir, "results", "submission_figures"); dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)

tidy <- function(x) {                       # round for readability; keep full precision in P values (3 significant digits)
  x <- copy(x)
  for (cn in names(x)) if (is.numeric(x[[cn]])) {
    if (grepl("(^|_)p$|pval|pvalue|min_.*_p$|^P$|mr_p|q_fdr|q_p", cn, ignore.case = TRUE)) x[, (cn) := signif(get(cn), 3)]
    else if (grepl("PP|pp|pip|r2|s_gwas|s_eqtl|score|raw|quantile", cn)) x[, (cn) := round(get(cn), 4)]
    else if (!all(x[[cn]] == round(x[[cn]]), na.rm = TRUE)) x[, (cn) := signif(get(cn), 4)] }
  x }
ord <- function(x) { k <- intersect(c("window", "locus", "ancestry"), names(x)); if ("PP.H4" %in% names(x)) setorderv(x, c(k, "PP.H4"), c(rep(1, length(k)), -1)) else if (length(k)) setorderv(x, k); x }

S <- list()
S$S1 <- ord(tidy(fread(file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"))))
S$S2 <- ord(tidy(fread(file.path(atl_dir, "Atlas_kidney_coloc_all.tsv.gz"))))
S$S3 <- tidy(fread(file.path(atl_dir, "Atlas_susie_coloc.tsv")))
S$S4 <- ord(tidy(fread(file.path(atl_dir, "Atlas_eqtlcatalogue_coloc.tsv.gz"))))
S$S5 <- ord(tidy(fread(file.path(atl_dir, "Atlas_JOB_pqtl_eqtl.tsv")))[, `:=`(MR_OR = round(exp(mr_beta), 3), MR_OR_lo = round(exp(mr_beta - 1.96 * mr_se), 3), MR_OR_hi = round(exp(mr_beta + 1.96 * mr_se), 3))])
S$S6 <- tidy(fread(file.path(atl_dir, "Aim1_ancestry_decomposition.tsv")))
ag <- fread(file.path(atl_dir, "alphagenome", "AlphaGenome_all_scores.tsv.gz"))
keep <- intersect(c("rsid", "locus", "variant_id", "gene_name", "variant_scorer", "output_type", "biosample_name", "biosample_type",
                    "Assay title", "junction_Start", "junction_End", "raw_score", "quantile_score", "risk_raw"), names(ag))
ag <- ag[, ..keep][, is_activity := grepl("active", variant_scorer, ignore.case = TRUE)]
fwrite(ag, file.path(sub_dir, "Supplemental_Table_S7_full_AlphaGenome.tsv.gz"), sep = "\t")      # complete scores as a separate file
grp <- intersect(c("rsid", "output_type", "is_activity"), names(ag))
S$S7 <- tidy(ag[order(-abs(raw_score))][, head(.SD, 50), by = grp])  # top 50 tracks per variant x output type
if (!file.exists(file.path(atl_dir, "Matched_testing_summary.tsv"))) stop("Run 25_matched_testing.R first (writes the S8 inputs to atl_dir).")
S$S8 <- tidy(fread(file.path(atl_dir, "Matched_testing_summary.tsv")))
S$S8_genes <- tidy(fread(file.path(atl_dir, "Matched_testing_gene_level.tsv")))
S$S8_iterations <- fread(file.path(atl_dir, "Matched_testing_iterations.tsv"))
S$S8_lineage_pairs <- fread(file.path(atl_dir, "Matched_testing_lineage_pairs.tsv"))
S$S8_strata <- fread(file.path(atl_dir, "Matched_testing_strata.tsv"))
cat(sprintf("%s: %d rows\n", names(S), sapply(S, nrow)), sep = "")

readme <- data.table(
  Sheet = c(names(S)),
  Title = c("Colocalization of IgAN GWAS signals with ImmuNexUT immune-cell eQTL (28 cell types; Japanese donors)",
            "Colocalization with kidney eQTL (NephQTL2 glomerulus and tubulointerstitium)",
            "SuSiE-based colocalization for prioritized immune-cell pairs",
            "Colocalization with eQTL Catalogue expression (ge) and splicing (LeafCutter) QTL: OneK1K, BLUEPRINT, Quach, Schmiedel, GTEx",
            "Colocalization and cis-Mendelian randomization with Japanese whole-blood eQTL and plasma pQTL (Japan COVID-19 Task Force; truncated at P < 0.05, approximate)",
            "Ancestry decomposition of the 30 IgAN signals (variance contribution on the log-odds scale, 2p(1-p)beta^2)",
            "AlphaGenome variant-effect scores for 24 non-HLA lead variants (top 50 tracks per variant and output type; complete scores in Supplemental_Table_S7_full_AlphaGenome.tsv.gz)",
            "Matched-testing sensitivity analysis, immune cells vs kidney: summary of 2,000 Monte Carlo draws per analysis",
            "Matched-testing sensitivity analysis: gene-level comparison restricted to genes with a strong eQTL in both resources",
            "Matched-testing sensitivity analysis: number of colocalized windows in each of the 2,000 draws",
            "Matched-testing sensitivity analysis: all 45 pairs of immune lineages",
            "Matched-testing sensitivity analysis: number of tests per window x GWAS ancestry stratum"),
  Notes = c("One row per window x ancestry (GWAS) x cell type x gene. PP.H3/PP.H4, posterior probabilities of distinct/shared causal variants (coloc.abf; p1 = p2 = 1e-4, p12 = 1e-5). top_snp_H4, variant with the highest posterior under H4 (GRCh37). dir_at_top, +1 if the allele raising expression raises IgAN risk.",
            "As S1. tissue: Glom, glomerulus; Tube, tubulointerstitium.",
            "coloc.susie with up to five signals per trait; s_gwas/s_eqtl, estimated LD-mismatch parameter; best_pair, credible-set leads of the best-supported pair.",
            "dataset, eQTL Catalogue accession; quant: ge, gene expression; leafcutter, splicing cluster (trait = intron coordinates, GRCh38). Splicing clusters tested only if min P < 1e-5.",
            "layer: eQTL (whole blood, n = 1,019) or pQTL (Olink Explore 3072, n = 1,384). MR, Wald ratio at the lead QTL variant (qtl_lead); single instrument, P equals the GWAS P at that variant; OR per unit of the normalized trait.",
            "f_eur/f_eas, risk-allele frequencies; b_*, log odds ratios; V_*, 2p(1-p)beta^2; d_freq/d_effect, exact split of V_EAS - V_EUR; q_fdr, Benjamini-Hochberg adjusted P for the effect difference.",
            "raw_score: signed log fold change for expression scorers (risk_raw re-expressed per risk allele); unsigned magnitude for splicing scorers; activity scorers (is_activity = TRUE) give predicted activity, not allelic effect. quantile_score, percentile relative to common variants for the same scorer and track.",
            "Restricted to the 23 windows tested in both ImmuNexUT and NephQTL2 (glomerulus and tubulointerstitium pooled). Cell-type matched: 2 of 28 immune-cell types drawn without replacement, equal probability. Lineage matched: 2 of 10 lineages, one cell type each. Test-count matched: within each window x GWAS ancestry, as many immune gene-cell-type tests drawn without replacement as kidney tests. A window counts if any sampled test reaches PP.H4 >= 0.8. empirical_P = (draws <= kidney count + 1)/(2,000 + 1). Seed 20260921 (script 25_matched_testing.R).",
            "Same gene, window and GWAS ancestry with eQTL P < 1e-5 in both resources; best PP.H4 across immune-cell types and across kidney compartments.",
            "cell_type_matched, test_count_matched and lineage_matched: windows with PP.H4 >= 0.8 in each draw.",
            "windows_colocalized: windows (of 23) with PP.H4 >= 0.8 in any cell type of either lineage.",
            "n_immune, ImmuNexUT tests; nk, NephQTL2 tests. Immune tests exceeded kidney tests in every stratum, so test-count matching was exact."))

wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", fgFill = "#E8F1FC", border = "bottom")
addWorksheet(wb, "README"); writeData(wb, "README", "Supplemental Tables S1-S8 for: Ancestry-Resolved Genetic Mapping Links IgA Nephropathy Risk to Myeloid APRIL, FCAR Splicing, and Circulating Factor H", startRow = 1)
writeData(wb, "README", readme, startRow = 3, headerStyle = hs); setColWidths(wb, "README", 1:3, c(8, 70, 110))
addStyle(wb, "README", createStyle(wrapText = TRUE, valign = "top"), rows = 4:(3 + nrow(readme)), cols = 2:3, gridExpand = TRUE)
for (s in names(S)) {
  addWorksheet(wb, s); writeData(wb, s, S[[s]], headerStyle = hs); freezePane(wb, s, firstRow = TRUE)
  setColWidths(wb, s, seq_along(S[[s]]), "auto") }
out <- file.path(sub_dir, "Supplemental_Tables_S1-S8.xlsx"); saveWorkbook(wb, out, overwrite = TRUE)
cat("Saved:", out, "\n"); system(paste("open", shQuote(path.expand(sub_dir))))
