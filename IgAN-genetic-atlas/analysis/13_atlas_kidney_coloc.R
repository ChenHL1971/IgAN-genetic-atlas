## ---- Step 13: kidney (NephQTL2) colocalization for all atlas windows (paste after Step 11b Part 1) ----
eq_all <- list.files(file.path(proj_dir, "data", "download_res"), pattern = "NephQTL2\\.txt\\.gz$", full.names = TRUE)
eqtl_files <- c(Glom = grep("Glom", eq_all, value = TRUE)[1], Tube = grep("Tub", eq_all, value = TRUE)[1]); print(basename(eqtl_files))
N_eqtl <- c(Glom = 240, Tube = 311)
kd_dir <- file.path(atl_dir, "coloc_kidney"); dir.create(kd_dir, showWarnings = FALSE)
# one pass per NephQTL2 file -> per-window caches (CHR = column 3, POS = column 4)
for (t in names(eqtl_files)) {
  outs <- file.path(kd_dir, sprintf("%s_eQTL_%s.tsv", win$wname, t))
  if (all(file.exists(outs))) { message("cached: ", t); next }
  cond <- paste(sprintf("($3==%d && $4>=%d && $4<=%d)", win$chr, win$start, win$end), collapse = " || ")
  message("Scanning NephQTL2 ", t, " (5-10 min) ...")
  x <- fread(cmd = sprintf("gzip -cd %s | awk -F'\\t' 'NR==1 || %s'", shQuote(path.expand(eqtl_files[[t]])), cond))
  setnames(x, make.names(names(x))); x[, gene := sub("\\..*$", "", gene)]
  for (i in seq_len(nrow(win))) fwrite(x[CHR == win$chr[i] & POS >= win$start[i] & POS <= win$end[i]], outs[i], sep = "\t")
  message("  ", t, ": ", nrow(x), " rows")
}
sym <- function(ids) { if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(setNames(ids, ids))
  s <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ids), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  out <- s[ids]; out[is.na(out)] <- ids[is.na(out)]; unname(out) }
kid_window <- function(w, min_snps = 50) {
  of <- file.path(kd_dir, paste0("coloc_", w, ".tsv")); if (file.exists(of)) return(fread(of))
  res <- list()
  for (t in names(eqtl_files)) {
    eq <- fread(file.path(kd_dir, sprintf("%s_eQTL_%s.tsv", w, t)))
    if (!nrow(eq)) next
    eq <- eq[!is.na(beta) & se > 0][, SNP := paste(CHR, POS, sep = ":")]
    for (a in c("EAS", "EUR")) {
      g <- get_gwas(w, a)
      m <- merge(eq, g[, .(SNP, gA1 = A1, gA2 = A2, BETA, varbeta, P)], by = "SNP")
      m[, dir := fifelse(gA1 == ALT & gA2 == REF, 1, fifelse(gA1 == REF & gA2 == ALT, -1, NA_real_))]
      m <- m[!is.na(dir)][, gwas_beta := BETA * dir][, MAF := pmin(AltFreq, 1 - AltFreq)][MAF >= 0.01]
      m <- unique(m, by = c("gene", "SNP"))
      for (gn in m[, .N, by = gene][N >= min_snps, gene]) {
        d <- m[gene == gn]
        invisible(capture.output(r <- suppressMessages(suppressWarnings(coloc.abf(
          list(beta = d$gwas_beta, varbeta = d$varbeta, snp = d$SNP, type = "cc", s = ncase[[a]] / (ncase[[a]] + nctrl[[a]])),
          list(beta = d$beta, varbeta = d$se^2, snp = d$SNP, type = "quant", N = N_eqtl[[t]], MAF = d$MAF))))))
        sm <- r$summary; top <- r$results[which.max(r$results$SNP.PP.H4), ]
        res[[length(res) + 1]] <- data.table(window = w, ancestry = a, tissue = t, gene_id = gn, nsnps = nrow(d),
          min_gwas_p = min(d$P), min_eqtl_p = min(d$p.value), PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]],
          top_snp_H4 = top$snp, top_SNP.PP.H4 = top$SNP.PP.H4)
      }
    }
  }
  res <- rbindlist(res); if (nrow(res)) { res[, gene := sym(gene_id)]; fwrite(res, of, sep = "\t") }
  message(sprintf("  kidney %s: %d tests", w, nrow(res))); res
}
kr <- rbindlist(lapply(win$wname, function(w) tryCatch(kid_window(w), error = function(e) { message("ERROR ", w, ": ", conditionMessage(e)); NULL })), fill = TRUE)
fwrite(kr, file.path(atl_dir, "Atlas_kidney_coloc_all.tsv.gz"), sep = "\t")
cat(sprintf("\nKidney tests: %d | PP.H4 >= 0.8: %d | 0.5-0.8: %d\n", nrow(kr), sum(kr$PP.H4 >= 0.8), sum(kr$PP.H4 >= 0.5 & kr$PP.H4 < 0.8)))
cat("\n[K1] kidney pairs with PP.H4 >= 0.5:\n")
print(kr[PP.H4 >= 0.5, .(window, ancestry, tissue, gene, nsnps, min_eqtl_p = signif(min_eqtl_p, 2), PP.H3 = round(PP.H3, 2),
      PP.H4 = round(PP.H4, 2), top_snp_H4, top_pp = round(top_SNP.PP.H4, 2))][order(-PP.H4)], nrows = 60)
imm <- fread(file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"))
cmp <- merge(kr[, .(kidney_maxPP = round(max(PP.H4), 2)), by = .(window, ancestry)],
             imm[, .(immune_maxPP = round(max(PP.H4), 2)), by = .(window, ancestry)], by = c("window", "ancestry"), all = TRUE)
cat("\n[K2] max PP.H4, kidney vs immune, per window:\n")
print(dcast(cmp, window ~ ancestry, value.var = c("kidney_maxPP", "immune_maxPP")), nrows = 40)
