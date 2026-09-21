# =============================================================================
# Paper 2 (atlas) - Step 26: ancestry-specific fine-mapping of the IgAN GWAS signal
#   at TNFSF12_13 (rs3803800), FCAR (rs1865097) and CFH (rs6677604).
#   SuSiE-RSS on each ancestry-specific GWAS with ancestry-matched 1000 Genomes LD
#   (503 EUR / 504 EAS), up to 5 signals; 95% credible sets.
# Uses the cached 1000G window files from Step 11b (results/atlas/kg/<window>.rds) and the
# GWAS window files from Step 10. Runs in a fresh R session (~2-5 min).
# Output: results/atlas/Finemap_GWAS_summary.tsv, Finemap_GWAS_credible_sets.tsv
# =============================================================================
proj_dir <- "~/Desktop/IgAN_Genetics_Paper2"
atl_dir  <- file.path(proj_dir, "results", "atlas")
win_dir  <- file.path(atl_dir, "windows")
kg_dir   <- file.path(atl_dir, "kg")
panel_f  <- file.path(proj_dir, "data", "1000G", "integrated_call_samples_v3.20130502.ALL.panel")
chain_f  <- file.path(proj_dir, "data", "hg38ToHg19.over.chain")
ncase <- c(EUR = 5556, EAS = 4590); nctrl <- c(EUR = 21178, EAS = 7573)

for (p in c("data.table", "susieR")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(susieR)
panel <- fread(panel_f, header = FALSE, fill = TRUE, select = 1:4, col.names = c("sample", "pop", "super_pop", "gender"))[sample != "sample"]

## lead variants (GRCh37). FCAR lead is given in GRCh38 and lifted here.
ch <- rtracklayer::import.chain(path.expand(chain_f))
l <- rtracklayer::liftOver(GenomicRanges::GRanges("chr19", IRanges::IRanges(54885805, 54885805)), ch)
fcar19 <- as.integer(BiocGenerics::start(unlist(l)))
loci <- data.table(window = c("TNFSF12_13", "FCAR", "CFH"), rsid = c("rs3803800", "rs1865097", "rs6677604"),
                   lead = c("17:7462969", paste0("19:", fcar19), "1:196686918"))
print(loci)

get_gwas <- function(w, a) {
  x <- fread(file.path(win_dir, sprintf("%s_GWAS_%s.tsv", w, a)))
  x <- x[!is.na(BETA) & SE > 0 & !is.na(P)][, `:=`(A1 = toupper(A1), A2 = toupper(A2), SNP = paste(CHR, BP_hg19, sep = ":"))]
  x <- x[!((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A") | (A1 == "C" & A2 == "G") | (A1 == "G" & A2 == "C"))]
  unique(x, by = "SNP")
}

finemap <- function(w, anc, lead, rsid, L = 5, flank = 250000) {
  k <- readRDS(file.path(kg_dir, paste0(w, ".rds")))
  smp <- intersect(panel[super_pop == anc, sample], colnames(k$D))
  lead_pos <- as.integer(sub(".*:", "", lead))
  g <- get_gwas(w, anc)
  m <- merge(g, k$meta, by = "SNP")[abs(POS19 - lead_pos) <= flank]
  m[, dir := fifelse(A1 == kALT & A2 == kREF, 1, fifelse(A1 == kREF & A2 == kALT, -1, NA_real_))]
  m <- m[!is.na(dir)][order(POS19)]
  X <- k$D[match(m$SNP, k$meta$SNP), smp, drop = FALSE]
  af <- rowMeans(X) / 2; ok <- af >= 0.01 & af <= 0.99
  m <- m[ok]; X <- X[ok, , drop = FALSE]
  R <- cor(t(X)); z <- m$BETA * m$dir / m$SE          # z per 1000G ALT allele, same orientation as R
  n <- ncase[[anc]] + nctrl[[anc]]
  s <- tryCatch(estimate_s_rss(z, R, n = n), error = function(e) NA_real_)
  fit <- susie_rss(z = z, R = R, n = n, L = L, coverage = 0.95, min_abs_corr = 0.5)
  pip <- setNames(fit$pip, m$SNP); cs <- fit$sets$cs
  lead_in <- lead %in% m$SNP
  top <- m[which.max(abs(z))]
  sm <- data.table(window = w, rsid = rsid, ancestry = anc, n_snps = nrow(m), s_mismatch = round(s, 3),
                   gwas_top = top$SNP, gwas_top_P = signif(top$P, 3), lead = lead, lead_in_data = lead_in,
                   lead_P = if (lead_in) signif(m[SNP == lead, P], 3) else NA_real_,
                   lead_PIP = if (lead_in) round(pip[[lead]], 3) else NA_real_,
                   n_cs = length(cs))
  csd <- if (length(cs)) rbindlist(lapply(seq_along(cs), function(i) {
    ix <- cs[[i]]; o <- ix[order(-pip[ix])]
    data.table(window = w, ancestry = anc, cs = i, cs_size = length(ix), purity_min_r = round(fit$sets$purity$min.abs.corr[i], 3),
               lead_in_cs = lead %in% names(pip)[ix], top_snp = names(pip)[o[1]], top_pip = round(pip[o[1]], 3),
               snps = paste(sprintf("%s(%.3f)", names(pip)[o], pip[o]), collapse = ";"))
  })) else data.table(window = w, ancestry = anc, cs = NA_integer_, cs_size = NA_integer_, purity_min_r = NA_real_,
                      lead_in_cs = NA, top_snp = NA_character_, top_pip = NA_real_, snps = NA_character_)
  list(summary = sm, cs = csd, pip = data.table(window = w, ancestry = anc, SNP = m$SNP, P = m$P, pip = round(fit$pip, 4)))
}

res <- lapply(seq_len(nrow(loci)), function(i) lapply(c("EAS", "EUR"), function(a) {
  message(sprintf("SuSiE-RSS %s %s", loci$window[i], a))
  tryCatch(finemap(loci$window[i], a, loci$lead[i], loci$rsid[i]), error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL })
}))
res <- unlist(res, recursive = FALSE); res <- res[!sapply(res, is.null)]
S  <- rbindlist(lapply(res, `[[`, "summary"), fill = TRUE)
CS <- rbindlist(lapply(res, `[[`, "cs"), fill = TRUE)
PIP <- rbindlist(lapply(res, `[[`, "pip"))

## overlap of 95% credible sets between ancestries (shared variants)
ov <- CS[!is.na(cs), .(snps = list(sub("\\(.*", "", strsplit(snps, ";")[[1]]))), by = .(window, ancestry, cs)]
ovl <- ov[, {
  e <- .SD[ancestry == "EAS"]; u <- .SD[ancestry == "EUR"]
  if (!nrow(e) || !nrow(u)) data.table(eas_cs = NA_integer_, eur_cs = NA_integer_, n_shared = NA_integer_) else
    CJ(i = seq_len(nrow(e)), j = seq_len(nrow(u)))[, .(eas_cs = e$cs[i], eur_cs = u$cs[j],
       n_shared = mapply(function(a, b) length(intersect(e$snps[[a]], u$snps[[b]])), i, j))]
}, by = window]

fwrite(S, file.path(atl_dir, "Finemap_GWAS_summary.tsv"), sep = "\t")
fwrite(CS, file.path(atl_dir, "Finemap_GWAS_credible_sets.tsv"), sep = "\t")
fwrite(PIP, file.path(atl_dir, "Finemap_GWAS_pip.tsv.gz"), sep = "\t")
cat("\n==== Summary ====\n"); print(S)
cat("\n==== Credible sets ====\n"); print(CS[, .(window, ancestry, cs, cs_size, purity_min_r, lead_in_cs, top_snp, top_pip)])
cat("\n==== EAS-EUR credible-set overlap ====\n"); print(ovl)

## ---- CFH sensitivity: single-signal model (L = 1). Its PIPs depend only on the z scores, not on the LD reference,
##      so it is robust to the summary-statistic/LD inconsistency seen with the European GWAS (s ~ 0.30; multi-signal fit
##      did not converge). The purity filter for credible sets still uses LD.
cfh1 <- lapply(c("EAS", "EUR"), function(a) finemap("CFH", a, loci[window == "CFH"]$lead, "rs6677604", L = 1))
S1 <- rbindlist(lapply(cfh1, `[[`, "summary"), fill = TRUE)[, model := "L=1"]
CS1 <- rbindlist(lapply(cfh1, `[[`, "cs"), fill = TRUE)[, model := "L=1"]
fwrite(S1, file.path(atl_dir, "Finemap_GWAS_CFH_L1_summary.tsv"), sep = "\t")
fwrite(CS1, file.path(atl_dir, "Finemap_GWAS_CFH_L1_credible_sets.tsv"), sep = "\t")
cat("\n==== CFH single-signal model ====\n"); print(S1); print(CS1[, .(ancestry, cs, cs_size, lead_in_cs, top_snp, top_pip)])
