## ---- Step 15: eQTL Catalogue — European single-cell eQTL replication, liver/kidney, and splicing QTL ----
## Paste after Step 11b Part 1. Remote tabix (hg38), no full downloads.
dsets <- fread(text = "
id,study,tissue,quant,n,url
QTD000606,OneK1K,B cell,ge,977,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000606/QTD000606.all.tsv.gz
QTD000607,OneK1K,memory B cell,ge,981,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000607/QTD000607.all.tsv.gz
QTD000609,OneK1K,monocyte,ge,959,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000609/QTD000609.all.tsv.gz
QTD000610,OneK1K,CD16+ monocyte,ge,930,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000610/QTD000610.all.tsv.gz
QTD000612,OneK1K,CD4+ T cell,ge,981,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000612/QTD000612.all.tsv.gz
QTD000620,OneK1K,NK cell,ge,981,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000620/QTD000620.all.tsv.gz
QTD000626,OneK1K,dendritic cell,ge,869,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000626/QTD000626.all.tsv.gz
QTD000629,OneK1K,plasmacytoid dendritic cell,ge,639,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000038/QTD000629/QTD000629.all.tsv.gz
QTD000021,BLUEPRINT,monocyte,ge,191,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000021/QTD000021.all.tsv.gz
QTD000026,BLUEPRINT,neutrophil,ge,196,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000026/QTD000026.all.tsv.gz
QTD000031,BLUEPRINT,CD4+ T cell,ge,167,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000031/QTD000031.all.tsv.gz
QTD000356,GTEx,blood,ge,670,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz
QTD000266,GTEx,liver,ge,208,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz
QTD000261,GTEx,kidney cortex,ge,73,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000261/QTD000261.all.tsv.gz
QTD000025,BLUEPRINT,monocyte,leafcutter,191,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000025/QTD000025.cc.tsv.gz
QTD000030,BLUEPRINT,neutrophil,leafcutter,196,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000030/QTD000030.cc.tsv.gz
QTD000035,BLUEPRINT,CD4+ T cell,leafcutter,167,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000035/QTD000035.cc.tsv.gz
QTD000413,Quach_2016,monocyte,leafcutter,200,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000024/QTD000413/QTD000413.cc.tsv.gz
QTD000478,Schmiedel_2018,B cell,leafcutter,91,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000478/QTD000478.cc.tsv.gz
QTD000508,Schmiedel_2018,monocyte,leafcutter,91,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000508/QTD000508.cc.tsv.gz
QTD000483,Schmiedel_2018,CD4+ T cell,leafcutter,88,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000483/QTD000483.cc.tsv.gz
QTD000360,GTEx,blood,leafcutter,670,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000360/QTD000360.cc.tsv.gz
QTD000270,GTEx,liver,leafcutter,208,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000270/QTD000270.cc.tsv.gz
QTD000265,GTEx,kidney cortex,leafcutter,73,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000265/QTD000265.cc.tsv.gz
QTD000225,GTEx,LCL,leafcutter,147,https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000225/QTD000225.cc.tsv.gz")
EC_COLS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se",
             "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
ec_dir <- file.path(atl_dir, "eqtlcat"); dir.create(ec_dir, showWarnings = FALSE)
w38 <- fread(file.path(atl_dir, "windows_hg38.tsv"))[name %in% win$wname]
w38[, chr := sub("^chr", "", chr)]

## fetch one dataset x one window (cached)
ec_fetch <- function(did, w) {
  f <- file.path(ec_dir, sprintf("%s_%s.rds", did, w)); if (file.exists(f)) return(readRDS(f))
  r <- w38[name == w]; u <- dsets[id == did, url]
  tf <- Rsamtools::TabixFile(u)
  ln <- tryCatch(Rsamtools::scanTabix(tf, param = GenomicRanges::GRanges(r$chr, IRanges::IRanges(r$start, r$end)))[[1]],
                 error = function(e) { message("  fetch failed ", did, " ", w, ": ", conditionMessage(e)); character(0) })
  x <- if (length(ln)) fread(text = ln, sep = "\t", header = FALSE) else data.table()
  if (nrow(x)) { if (ncol(x) == length(EC_COLS)) setnames(x, EC_COLS) else stop("unexpected column count ", ncol(x), " in ", did)
    x <- x[type == "SNP" | (nchar(ref) == 1 & nchar(alt) == 1), .(molecular_trait_id, gene_id, position, ref, alt, maf, pvalue, beta, se)] }
  saveRDS(x, f); x
}
## coloc of one dataset x window x ancestry
ec_coloc <- function(did, w, a, min_snps = 50) {
  x <- ec_fetch(did, w); if (!nrow(x)) return(NULL)
  q <- dsets[id == did]
  if (q$quant == "leafcutter") { keep <- x[, .(mp = min(pvalue)), by = molecular_trait_id][mp < 1e-5, molecular_trait_id]; x <- x[molecular_trait_id %in% keep] }
  if (!nrow(x)) return(NULL)
  up <- unique(x[, .(position)])[, pos19 := lift19(paste0("chr", w38[name == w, chr]), position)]
  x <- merge(x, up[!is.na(pos19)], by = "position")[, SNP := paste(w38[name == w, chr], pos19, sep = ":")]
  g <- get_gwas(w, a)
  m <- merge(x, g[, .(SNP, gA1 = A1, gA2 = A2, BETA, varbeta, P)], by = "SNP")
  m[, dir := fifelse(gA1 == alt & gA2 == ref, 1, fifelse(gA1 == ref & gA2 == alt, -1, NA_real_))]
  m <- m[!is.na(dir) & maf >= 0.01 & se > 0][, gwas_beta := BETA * dir]
  m <- unique(m, by = c("molecular_trait_id", "SNP"))
  tr <- m[, .N, by = .(molecular_trait_id, gene_id)][N >= min_snps]
  rbindlist(lapply(seq_len(nrow(tr)), function(i) {
    d <- m[molecular_trait_id == tr$molecular_trait_id[i]]
    invisible(capture.output(r <- suppressMessages(suppressWarnings(coloc.abf(
      list(beta = d$gwas_beta, varbeta = d$varbeta, snp = d$SNP, type = "cc", s = ncase[[a]] / (ncase[[a]] + nctrl[[a]])),
      list(beta = d$beta, varbeta = d$se^2, snp = d$SNP, type = "quant", N = q$n, MAF = d$maf))))))
    sm <- r$summary; top <- r$results[which.max(r$results$SNP.PP.H4), ]; ti <- match(top$snp, d$SNP)
    data.table(dataset = did, study = q$study, tissue = q$tissue, quant = q$quant, window = w, ancestry = a,
               trait = tr$molecular_trait_id[i], gene_id = tr$gene_id[i], nsnps = nrow(d), min_qtl_p = min(d$pvalue),
               PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]], top_snp_H4 = top$snp, top_pp = top$SNP.PP.H4,
               dir_at_top = sign(d$beta[ti]) * sign(d$gwas_beta[ti]))
  }))
}
## run: all datasets x all windows, European GWAS first (ancestry-matched), then East Asian
t0 <- Sys.time(); ecr <- list()
for (k in seq_len(nrow(dsets))) for (w in w38$name) {
  did <- dsets$id[k]; message(sprintf("%s (%s %s %s) x %s", did, dsets$study[k], dsets$tissue[k], dsets$quant[k], w))
  for (a in c("EUR", "EAS")) ecr[[length(ecr) + 1]] <- tryCatch(ec_coloc(did, w, a), error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL })
}
ecr <- rbindlist(ecr, fill = TRUE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
ecr[, gene := unname(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))]
ecr[is.na(gene), gene := gene_id]
fwrite(ecr, file.path(atl_dir, "Atlas_eqtlcatalogue_coloc.tsv.gz"), sep = "\t")
message(sprintf("done: %d tests in %.0f min", nrow(ecr), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

targets <- c("TNFSF13", "LYN", "CD28", "CARD9", "ITGAM", "ITGAX", "FCAR", "CFH", "CFHR1", "FCRL5", "FCRL4", "IRF8", "TNFSF4", "OSM", "LIF")
cat("\n[E1] Target genes, best result per gene x dataset (PP.H4):\n")
tt <- ecr[gene %in% targets][order(-PP.H4)][, .SD[1], by = .(gene, study, tissue, quant, ancestry)][, ds := paste(study, tissue, quant, sep = "|")]
print(dcast(tt, gene + ancestry ~ ds, value.var = "PP.H4", fun.aggregate = function(x) round(x[1], 2)), nrows = 60)
cat("\n[E2] All pairs with PP.H4 >= 0.8:\n")
print(ecr[PP.H4 >= 0.8, .(window, ancestry, study, tissue, quant, gene, nsnps, min_qtl_p = signif(min_qtl_p, 2), PP.H4 = round(PP.H4, 2),
      top_snp_H4, top_pp = round(top_pp, 2), dir_at_top)][order(window, -PP.H4)], nrows = 150)
