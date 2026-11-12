## ---- Step 27: sensitivity of coloc.abf to truncating QTL data at P < 0.05 ----
## The Japan COVID-19 Task Force whole-blood eQTL and plasma pQTL files keep only variant-trait pairs with
## P < 0.05 (or PIP > 0.001). This script asks how much that truncation changes PP.H4, by applying the
## same truncation to COMPLETE QTL data and re-running coloc.abf:
##   A. ImmuNexUT (complete nominal eQTL): every window x ancestry x cell x gene with published PP.H4 >= 0.5,
##      plus a random sample of 2,000 tests with PP.H4 < 0.5 (to check false positives).
##   B. eQTL Catalogue (complete): GTEx whole blood TNFSF13 and CFH, GTEx liver CFHR3/CFHR4, and the FCAR
##      junction 19:54885525:54885805 in Quach monocytes and GTEx whole blood.
## Truncated input = shared variants with QTL P < 0.05; tests with < 30 such variants are reported as NA,
## exactly as in Step 16.
##
## HOW TO RUN: paste after Step 11b Part 1 (defines get_gwas, kg_window, kg_freq, lift19, co_dir, ncase, nctrl,
## N_imx) and after the function definitions of Step 15 (dsets, ec_fetch, w38). Cached inputs are reused, so
## no new downloads are needed except for any eQTL Catalogue window not already cached.
## Outputs: results/atlas/truncation/truncation_immunexut.tsv and truncation_eqtlcatalogue.tsv (summarized in Supplementary Table S12).

set.seed(20261111)
tr_dir <- file.path(atl_dir, "truncation"); dir.create(tr_dir, showWarnings = FALSE)

run_coloc <- function(gb, gvb, snp, a, qb, qvb, N, maf) {
  invisible(capture.output(r <- suppressMessages(suppressWarnings(coloc.abf(
    list(beta = gb, varbeta = gvb, snp = snp, type = "cc", s = ncase[[a]] / (ncase[[a]] + nctrl[[a]])),
    list(beta = qb, varbeta = qvb, snp = snp, type = "quant", N = N, MAF = maf))))))
  top <- r$results[which.max(r$results$SNP.PP.H4), ]
  list(PP.H3 = r$summary[["PP.H3.abf"]], PP.H4 = r$summary[["PP.H4.abf"]], top = top$snp)
}

## full vs truncated for one prepared data.table d (columns: SNP, gwas_beta, varbeta, qb, qse, qp, MAF)
compare_one <- function(d, a, N) {
  f <- run_coloc(d$gwas_beta, d$varbeta, d$SNP, a, d$qb, d$qse^2, N, d$MAF)
  t <- d[qp < 0.05]
  if (nrow(t) >= 30) {
    tt <- run_coloc(t$gwas_beta, t$varbeta, t$SNP, a, t$qb, t$qse^2, N, t$MAF)
  } else tt <- list(PP.H3 = NA_real_, PP.H4 = NA_real_, top = NA_character_)
  data.table(nsnps_full = nrow(d), nsnps_trunc = nrow(t),
             PP.H4_full = f$PP.H4, PP.H4_trunc = tt$PP.H4, PP.H3_full = f$PP.H3, PP.H3_trunc = tt$PP.H3,
             top_full = f$top, top_trunc = tt$top)
}

## ---------- A. ImmuNexUT ----------
imx_prepared <- function(w) {   # same filtering as imx_window() in Step 11b
  e <- fread(file.path(co_dir, paste0("imx_", w, "_hg38.tsv")))
  e <- e[!is.na(slope.ALT.) & nominal_P_value > 0 & nominal_P_value < 1 & nchar(REF) == 1 & nchar(ALT) == 1]
  e[, z := qnorm(nominal_P_value / 2, lower.tail = FALSE) * sign(slope.ALT.)]; e <- e[abs(z) > 0][, se := abs(slope.ALT. / z)]
  up <- unique(e[, .(Variant_CHR, Variant_position_start)])[, pos19 := lift19(Variant_CHR, Variant_position_start)]
  e <- merge(e, up[!is.na(pos19)], by = c("Variant_CHR", "Variant_position_start"))
  e[, SNP := paste(sub("^chr", "", Variant_CHR), pos19, sep = ":")]
  e <- merge(e, kg_freq(kg_window(w)), by = "SNP")
  e[, kdir := fifelse(kREF == REF & kALT == ALT, 1, fifelse(kREF == ALT & kALT == REF, -1, NA_real_))]
  e <- e[!is.na(kdir)][, af_alt := fifelse(kdir == 1, AF_JPT, 1 - AF_JPT)][, MAF := pmin(af_alt, 1 - af_alt)][MAF >= 0.01]
  unique(e, by = c("cell", "Gene_id", "SNP"))
}
pub <- rbindlist(lapply(win$wname, function(w) { f <- file.path(co_dir, paste0("coloc_", w, ".tsv")); if (file.exists(f)) fread(f) }))
sel <- rbind(pub[PP.H4 >= 0.5], pub[PP.H4 < 0.5][sample(.N, min(.N, 2000))])
resA <- list()
for (w in unique(sel$window)) {
  e <- imx_prepared(w)
  for (a in c("EAS", "EUR")) {
    s <- sel[window == w & ancestry == a]; if (!nrow(s)) next
    g <- get_gwas(w, a)
    m <- merge(e, g[, .(SNP, gA1 = A1, gA2 = A2, BETA, varbeta)], by = "SNP")
    m[, gdir := fifelse(gA1 == ALT & gA2 == REF, 1, fifelse(gA1 == REF & gA2 == ALT, -1, NA_real_))]
    m <- m[!is.na(gdir)][, gwas_beta := BETA * gdir]
    for (i in seq_len(nrow(s))) {
      d <- m[cell == s$cell[i] & Gene_id == s$gene_id[i], .(SNP, gwas_beta, varbeta, qb = slope.ALT., qse = se, qp = nominal_P_value, MAF)]
      if (nrow(d) < 50) next
      resA[[length(resA) + 1]] <- cbind(data.table(source = "ImmuNexUT", window = w, ancestry = a, cell = s$cell[i], gene = s$gene[i],
                                                   PP.H4_published = s$PP.H4[i]), compare_one(d, a, N_imx))
    }
    message(sprintf("A  %s %s: %d tests", w, a, nrow(s)))
  }
}
resA <- rbindlist(resA)

## ---------- B. eQTL Catalogue targets ----------
targets <- data.table(
  did   = c("QTD000356", "QTD000356", "QTD000266", "QTD000266", "QTD000413", "QTD000360"),
  win   = c("TNFSF12_13", "CFH", "CFH", "CFH", "FCAR", "FCAR"),
  trait = c("ENSG00000161955", "ENSG00000000971", "ENSG00000116785", "ENSG00000134365",
            "19:54885525:54885805", "19:54885525:54885805"),
  label = c("GTEx blood TNFSF13", "GTEx blood CFH", "GTEx liver CFHR3", "GTEx liver CFHR4",
            "Quach monocyte FCAR junction", "GTEx blood FCAR junction"))
resB <- list()
for (k in seq_len(nrow(targets))) {
  x <- ec_fetch(targets$did[k], targets$win[k]); if (!nrow(x)) next
  x <- x[grepl(targets$trait[k], molecular_trait_id, fixed = TRUE)]
  up <- unique(x[, .(position)])[, pos19 := lift19(paste0("chr", w38[name == targets$win[k], chr]), position)]
  x <- merge(x, up[!is.na(pos19)], by = "position")[, SNP := paste(w38[name == targets$win[k], chr], pos19, sep = ":")]
  for (a in c("EAS", "EUR")) {
    g <- get_gwas(targets$win[k], a)
    m <- merge(x, g[, .(SNP, gA1 = A1, gA2 = A2, BETA, varbeta)], by = "SNP")
    m[, dir := fifelse(gA1 == alt & gA2 == ref, 1, fifelse(gA1 == ref & gA2 == alt, -1, NA_real_))]
    m <- unique(m[!is.na(dir) & maf >= 0.01 & se > 0][, gwas_beta := BETA * dir], by = "SNP")
    for (tid in unique(m$molecular_trait_id)) {
      d <- m[molecular_trait_id == tid, .(SNP, gwas_beta, varbeta, qb = beta, qse = se, qp = pvalue, MAF = maf)]
      if (nrow(d) < 50) next
      resB[[length(resB) + 1]] <- cbind(data.table(source = "eQTL Catalogue", label = targets$label[k], trait = tid, ancestry = a),
                                        compare_one(d, a, dsets[id == targets$did[k], n]))
    }
  }
}
resB <- rbindlist(resB)

## ---------- summaries ----------
resA[, `:=`(strong_full = PP.H4_full >= 0.8, strong_trunc = PP.H4_trunc >= 0.8, same_top = top_full == top_trunc)]
cat("\nA. ImmuNexUT, full vs truncated at P < 0.05\n")
print(resA[, .(tests = .N, trunc_NA = sum(is.na(PP.H4_trunc)),
               spearman = round(cor(PP.H4_full, PP.H4_trunc, method = "spearman", use = "complete.obs"), 3),
               median_abs_diff = round(median(abs(PP.H4_full - PP.H4_trunc), na.rm = TRUE), 3),
               strong_kept = sum(strong_full & strong_trunc, na.rm = TRUE), strong_lost = sum(strong_full & !strong_trunc, na.rm = TRUE),
               strong_gained = sum(!strong_full & strong_trunc, na.rm = TRUE), same_top_snp = mean(same_top, na.rm = TRUE))])
cat("\nA. Prioritized genes\n")
print(resA[gene %in% c("TNFSF13", "CARD9", "ITGAM", "CD28", "LYN", "FCRL5", "IRF8", "TNFSF4", "CFH"),
           .SD[which.max(PP.H4_full)], by = .(gene, ancestry)][, .(gene, ancestry, cell, PP.H4_full = round(PP.H4_full, 3),
           PP.H4_trunc = round(PP.H4_trunc, 3), nsnps_full, nsnps_trunc, same_top)])
cat("\nB. eQTL Catalogue targets\n")
print(resB[, .(label, trait, ancestry, nsnps_full, nsnps_trunc, PP.H4_full = round(PP.H4_full, 3), PP.H4_trunc = round(PP.H4_trunc, 3),
               same_top = top_full == top_trunc)])
fwrite(resA, file.path(tr_dir, "truncation_immunexut.tsv"), sep = "\t")
fwrite(resB, file.path(tr_dir, "truncation_eqtlcatalogue.tsv"), sep = "\t")
