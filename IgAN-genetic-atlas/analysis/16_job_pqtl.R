## ---- Step 16: Japanese plasma pQTL (Olink 3072, n=1,384) + whole-blood eQTL (n=1,019) from JOB / Wang 2024 ----
## Paste after Step 11b Part 1. Files in data/JOB/: pqtl_sumstats_..._taskforce.tsv.gz, eqtl_sumstats_..._taskforce.tsv.gz
## Note: files contain only pairs with P < 0.05 or PIP > 0.001 (truncated) -> coloc is approximate; PIPs used for fine-map overlap.
job_dir <- file.path(proj_dir, "data", "JOB"); jb_out <- file.path(atl_dir, "job"); dir.create(jb_out, showWarnings = FALSE)
job_files <- c(pQTL = list.files(job_dir, "^pqtl_sumstats.*\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)[1],
               eQTL = list.files(job_dir, "^eqtl_sumstats.*\\.tsv\\.gz$", full.names = TRUE, recursive = TRUE)[1]); print(basename(job_files))
N_job <- c(pQTL = 1384, eQTL = 1019)
## 1. one pass per file: keep rows whose hg19 variant lies in an atlas window (awk on column 1 "chr:pos:ref:alt")
for (k in names(job_files)) {
  of <- file.path(jb_out, paste0(k, "_windows.tsv.gz")); if (file.exists(of)) { message("cached: ", k); next }
  cond <- paste(sprintf("(c==\"%d\" && p>=%d && p<=%d)", win$chr, win$start, win$end), collapse = " || ")
  awk <- sprintf("NR==1 {print; next} {split($1, a, \":\"); c=a[1]; sub(/^chr/, \"\", c); p=a[2]+0; if (%s) print}", cond)
  message("Scanning ", k, " (about 5-15 min) ...")
  x <- fread(cmd = sprintf("gzip -cd %s | awk -F'\\t' '%s'", shQuote(path.expand(job_files[[k]])), awk))
  fwrite(x, of, sep = "\t"); message("  ", k, ": ", nrow(x), " rows in windows")
}
jb <- rbindlist(lapply(names(job_files), function(k) fread(file.path(jb_out, paste0(k, "_windows.tsv.gz")))[, layer := k]), fill = TRUE)
jb[, c("chr19", "pos19", "ref", "alt") := tstrsplit(variant_id_hg19, ":", fixed = TRUE)]
jb[, `:=`(chr19 = sub("^chr", "", chr19), pos19 = as.integer(pos19))]
jb <- jb[nchar(ref) == 1 & nchar(alt) == 1 & slope_se > 0][, SNP := paste(chr19, pos19, sep = ":")]
jb[, window := NA_character_]
for (i in seq_len(nrow(win))) jb[chr19 == as.character(win$chr[i]) & pos19 >= win$start[i] & pos19 <= win$end[i], window := win$wname[i]]
jb <- jb[!is.na(window)]
cat("\nProteins (pQTL) available per window:\n"); print(jb[layer == "pQTL", .(proteins = paste(sort(unique(gene_name)), collapse = ", ")), by = window], nrows = 30)

## 2. coloc (approximate, truncated summary statistics) per window x layer x gene x ancestry
jres <- list()
for (w in unique(jb$window)) for (k in c("pQTL", "eQTL")) for (a in c("EAS", "EUR")) {
  g <- get_gwas(w, a)
  m <- merge(jb[window == w & layer == k], g[, .(SNP, gA1 = A1, gA2 = A2, BETA, SE, varbeta, P)], by = "SNP")
  m[, dir := fifelse(gA1 == alt & gA2 == ref, 1, fifelse(gA1 == ref & gA2 == alt, -1, NA_real_))]
  m <- unique(m[!is.na(dir) & maf >= 0.01][, gwas_beta := BETA * dir], by = c("gene_id", "SNP"))
  for (gn in m[, .N, by = gene_name][N >= 30, gene_name]) {
    d <- m[gene_name == gn]
    invisible(capture.output(r <- suppressMessages(suppressWarnings(coloc.abf(
      list(beta = d$gwas_beta, varbeta = d$varbeta, snp = d$SNP, type = "cc", s = ncase[[a]] / (ncase[[a]] + nctrl[[a]])),
      list(beta = d$slope, varbeta = d$slope_se^2, snp = d$SNP, type = "quant", N = N_job[[k]], MAF = d$maf))))))
    sm <- r$summary; top <- r$results[which.max(r$results$SNP.PP.H4), ]; ti <- match(top$snp, d$SNP)
    lq <- d[which.min(pval_nominal)]; fm <- d[which.max(pip_susie)]
    jres[[length(jres) + 1]] <- data.table(window = w, layer = k, ancestry = a, gene = gn, nsnps = nrow(d),
      min_qtl_p = min(d$pval_nominal), PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]], top_snp_H4 = top$snp, top_pp = top$SNP.PP.H4,
      dir_at_top = sign(d$slope[ti]) * sign(d$gwas_beta[ti]),
      qtl_lead = lq$SNP, qtl_lead_rsid = lq$rsid, top_pip_snp = fm$SNP, top_pip = fm$pip_susie, gwas_p_at_top_pip = fm$P,
      ## cis-MR (Wald ratio) at the lead QTL variant: log-odds IgAN per 1-unit (NPX / normalised expression) increase
      mr_beta = lq$gwas_beta / lq$slope, mr_se = lq$SE / abs(lq$slope), mr_p = 2 * pnorm(-abs((lq$gwas_beta / lq$slope) / (lq$SE / abs(lq$slope)))))
  }
}
jres <- rbindlist(jres); fwrite(jres, file.path(atl_dir, "Atlas_JOB_pqtl_eqtl.tsv"), sep = "\t")
targets <- c("TNFSF13", "LYN", "CD28", "CARD9", "ITGAM", "ITGAX", "FCAR", "CFH", "CFHR1", "FCRL5", "FCRL4", "IRF8", "TNFSF4", "OSM", "LIF", "TNFRSF13B", "TNFSF8", "CFHR3", "CFHR4", "CXCL1", "RNASEH2C")
cat("\n[J1] Target genes (both layers, both ancestries):\n")
print(jres[gene %in% targets, .(layer, gene, ancestry, nsnps, min_qtl_p = signif(min_qtl_p, 2), PP.H3 = round(PP.H3, 2), PP.H4 = round(PP.H4, 2),
      top_snp_H4, qtl_lead, top_pip_snp, top_pip = round(top_pip, 2), dir_at_top, MR_OR = round(exp(mr_beta), 2), MR_P = signif(mr_p, 2))][order(layer, gene, ancestry)], nrows = 80)
cat("\n[J2] Any protein (pQTL) with PP.H4 >= 0.5:\n")
print(jres[layer == "pQTL" & PP.H4 >= 0.5, .(window, ancestry, gene, min_qtl_p = signif(min_qtl_p, 2), PP.H4 = round(PP.H4, 2), top_snp_H4,
      top_pip_snp, top_pip = round(top_pip, 2), dir_at_top, MR_OR = round(exp(mr_beta), 2), MR_P = signif(mr_p, 2))][order(-PP.H4)], nrows = 60)
