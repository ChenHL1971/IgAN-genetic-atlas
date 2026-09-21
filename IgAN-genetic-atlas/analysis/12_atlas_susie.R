## ---- Step 12: SuSiE-coloc for the strong atlas hits (paste AFTER Parts 1-2 of Step 11b) ----
for (p in c("susieR")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(susieR)
eas_s <- panel[super_pop == "EAS", sample]; eur_s <- panel[super_pop == "EUR", sample]
prep_pair <- function(w, cell_, gene_, anc) {
  e <- fread(file.path(co_dir, paste0("imx_", w, "_hg38.tsv")))[cell == cell_ & Gene_name == gene_]
  e <- e[!is.na(slope.ALT.) & nominal_P_value > 0 & nominal_P_value < 1 & nchar(REF) == 1 & nchar(ALT) == 1]
  e[, z := qnorm(nominal_P_value / 2, lower.tail = FALSE) * sign(slope.ALT.)]; e <- e[abs(z) > 0][, se := abs(slope.ALT. / z)]
  up <- unique(e[, .(Variant_CHR, Variant_position_start)])[, pos19 := lift19(Variant_CHR, Variant_position_start)]
  e <- merge(e, up[!is.na(pos19)], by = c("Variant_CHR", "Variant_position_start"))
  e[, SNP := paste(sub("^chr", "", Variant_CHR), pos19, sep = ":")]
  k <- kg_window(w); e <- merge(e, kg_freq(k), by = "SNP")
  e[, kdir := fifelse(kREF == REF & kALT == ALT, 1, fifelse(kREF == ALT & kALT == REF, -1, NA_real_))]
  e <- e[!is.na(kdir)][, af_alt := fifelse(kdir == 1, AF_JPT, 1 - AF_JPT)][, MAF := pmin(af_alt, 1 - af_alt)][MAF >= 0.01]
  e <- unique(e, by = "SNP")
  g <- get_gwas(w, anc)
  m <- merge(e, g[, .(SNP, gA1 = A1, gA2 = A2, BETA, varbeta, P)], by = "SNP")
  m[, gdir := fifelse(gA1 == ALT & gA2 == REF, 1, fifelse(gA1 == REF & gA2 == ALT, -1, NA_real_))]
  m <- m[!is.na(gdir)]
  # both effects re-expressed per 1000G ALT allele (same orientation as the LD matrix)
  m[, `:=`(b_g = BETA * gdir * kdir, b_e = slope.ALT. * kdir)]
  m[order(pos19)]
}
ld_of <- function(w, snps, samples) {
  k <- kg_window(w); i <- match(snps, k$meta$SNP)
  X <- k$D[i, intersect(samples, colnames(k$D)), drop = FALSE]
  ok <- apply(X, 1, sd) > 0
  R <- suppressWarnings(cor(t(X[ok, , drop = FALSE]))); dimnames(R) <- list(snps[ok], snps[ok]); R
}
susie_pair <- function(w, cell_, gene_, anc, L = 5) {
  d <- prep_pair(w, cell_, gene_, anc)
  Rg <- ld_of(w, d$SNP, if (anc == "EUR") eur_s else eas_s); Re <- ld_of(w, d$SNP, eas_s)
  keep <- intersect(rownames(Rg), rownames(Re)); d <- d[SNP %in% keep]; Rg <- Rg[d$SNP, d$SNP]; Re <- Re[d$SNP, d$SNP]
  zg <- d$b_g / sqrt(d$varbeta); ze <- d$b_e / d$se
  s_g <- tryCatch(estimate_s_rss(zg, Rg, n = ncase[[anc]] + nctrl[[anc]]), error = function(e) NA)
  s_e <- tryCatch(estimate_s_rss(ze, Re, n = N_imx), error = function(e) NA)
  D1 <- list(beta = setNames(d$b_g, d$SNP), varbeta = setNames(d$varbeta, d$SNP), snp = d$SNP, position = d$pos19,
             type = "cc", s = ncase[[anc]] / (ncase[[anc]] + nctrl[[anc]]), N = ncase[[anc]] + nctrl[[anc]], LD = Rg)
  D2 <- list(beta = setNames(d$b_e, d$SNP), varbeta = setNames(d$se^2, d$SNP), snp = d$SNP, position = d$pos19,
             type = "quant", N = N_imx, MAF = d$MAF, LD = Re)
  S1 <- tryCatch(suppressMessages(runsusie(D1, L = L)), error = function(e) NULL)
  S2 <- tryCatch(suppressMessages(runsusie(D2, L = L)), error = function(e) NULL)
  cs_top <- function(S) if (is.null(S) || is.null(S$sets$cs)) NA_character_ else
    paste(sapply(S$sets$cs, function(ix) names(S$pip)[ix][which.max(S$pip[ix])]), collapse = ";")
  base <- data.table(window = w, cell = cell_, gene = gene_, ancestry = anc, nsnps = nrow(d),
                     s_gwas = round(s_g, 3), s_eqtl = round(s_e, 3), gwas_cs_tops = cs_top(S1), eqtl_cs_tops = cs_top(S2))
  if (is.null(S1$sets$cs) || is.null(S2$sets$cs)) return(base[, `:=`(best_PP.H4 = NA_real_, best_pair = "no credible set in one trait")])
  sm <- as.data.table(coloc.susie(S1, S2)$summary)
  b <- sm[which.max(PP.H4.abf)]
  base[, `:=`(best_PP.H4 = round(b$PP.H4.abf, 3), best_PP.H3 = round(b$PP.H3.abf, 3), best_pair = paste(b$hit1, b$hit2, sep = " / "),
              n_cs_pairs = nrow(sm), n_pairs_H4gt0.8 = sum(sm$PP.H4.abf > 0.8))]
}
targets <- fread(text = "
window,cell,gene,ancestry
TNFSF12_13,CL_Mono,TNFSF13,EAS
LYN,Plasmablast,LYN,EAS
LYN,Plasmablast,LYN,EUR
LYN,SM_B,LYN,EUR
CD28,Naive_CD4,CD28,EAS
CARD9,NC_Mono,CARD9,EUR
CARD9,CL_Mono,CARD9,EUR
FCRL3,DN_B,FCRL5,EAS
FCRL3,SM_B,FCRL4,EAS
IRF8,pDC,IRF8,EUR
TNFSF4_18,Naive_B,TNFSF4,EUR
OVOL1_RELA,CL_Mono,AP5B1,EAS
ITGAM_ITGAX,Neu,ITGAX,EAS
ITGAM_ITGAX,CL_Mono,ITGAM,EAS")
sres <- rbindlist(lapply(seq_len(nrow(targets)), function(i) {
  t <- targets[i]; message(sprintf("SuSiE %s %s %s %s", t$window, t$cell, t$gene, t$ancestry))
  tryCatch(susie_pair(t$window, t$cell, t$gene, t$ancestry), error = function(e) data.table(window = t$window, gene = t$gene, best_pair = paste("ERROR:", conditionMessage(e))))
}), fill = TRUE)
fwrite(sres, file.path(atl_dir, "Atlas_susie_coloc.tsv"), sep = "\t")
print(sres[, .(window, cell, gene, ancestry, nsnps, s_gwas, s_eqtl, best_PP.H4, best_PP.H3, best_pair, n_pairs_H4gt0.8)])
