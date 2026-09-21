# =============================================================================
# Paper 2 (atlas) - Step 10: Aim 1 ancestry decomposition + locus windows for Aims 2-3
#   Part A  30 independent IgAN signals (Kiryluk 2023 Table 1)
#   Part B  one pass through each GWAS -> +/-500 kb windows cached per locus
#   Part C  Aim 1: EUR vs EAS effect, Cochran Q, frequency-vs-effect decomposition
#   Part D  ancestry-specific leads (+/-250 kb) looked up in the other ancestry
#   Part E  windows in hg38 (for ImmuNexUT / Japan Omics Browser, GRCh38)
# All positions hg19 unless stated.
# =============================================================================
proj_dir <- "~/Desktop/IgAN_Genetics_Paper2"
data_dir <- file.path(proj_dir, "data")
atl_dir  <- file.path(proj_dir, "results", "atlas"); dir.create(atl_dir, showWarnings = FALSE, recursive = TRUE)
win_dir  <- file.path(atl_dir, "windows");          dir.create(win_dir, showWarnings = FALSE)
gwas_files <- c(EUR = file.path(data_dir, "IgAN_Combined_metaanalysis_European_only.txt"),
                EAS = file.path(data_dir, "IgAN_Combined_metaanalysis_Asian_only.txt"))
FLANK <- 500000L

for (p in c("data.table", "ggplot2", "ggrepel", "R.utils")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(ggplot2)

## ---- Part A: the 30 signals (risk allele; control frequencies as published) ----
sig <- fread(text = "
chr,pos,rsid,locus,risk,f_eur,f_eas
1,157542162,rs849815,FCRL3,A,0.66,0.50
1,173146357,rs4916312,TNFSF4/18,A,0.35,0.07
1,196686918,rs6677604,CFH,G,0.80,0.93
1,196603302,rs12029571,CFH (2),A,0.22,0.35
2,61092678,rs842638,REL,T,0.44,0.15
2,204584759,rs3769684,CD28,T,0.95,0.49
4,74725320,rs6828610,PF4V1,G,0.16,0.27
6,249571,rs12201499,IRF4/DUSP22,C,0.12,0.28
6,7214676,rs12530084,LY86,C,0.77,0.51
6,32389305,rs9268557,HLA-DRA,C,0.51,0.57
6,32667829,rs9275355,HLA-DQB/DQA,C,0.23,0.33
6,32599999,rs9272105,HLA-DQA,A,0.60,0.45
6,32681631,rs9275596,HLA-DQB/DQA (2),T,0.66,0.81
6,33074288,rs3128927,HLA-DPA/DPB,C,0.73,0.83
8,6808722,rs2075836,DEFA1/4,T,0.31,0.30
8,56852496,rs75413466,LYN,A,0.02,0.06
8,124765474,rs34354351,ANXA3,T,0.17,0.32
9,117643362,rs13300483,TNFSF8/15,T,0.24,0.31
9,139266496,rs4077515,CARD9,T,0.41,0.29
10,65363048,rs57917667,REEP3,G,0.02,0.19
10,81043743,rs1108618,ZMIZ1,A,0.60,0.49
11,65555524,rs10896045,OVOL1/RELA,A,0.30,0.48
11,128487069,rs7121743,ETS1,C,0.16,0.47
14,107222014,rs751081288,IGH,A,0.43,0.53
16,31357760,rs11150612,ITGAM/ITGAX,A,0.64,0.27
16,86017715,rs1879210,IRF8,T,0.64,0.86
17,7462969,rs3803800,TNFSF12/13,A,0.21,0.32
17,16851450,rs57382045,TNFRSF13B,A,0.11,0.33
19,55397217,rs1865097,FCAR,A,0.30,0.38
22,30512478,rs4823074,LIF/OSM,G,0.54,0.67")
sig[, hla := chr == 6 & pos > 28.4e6 & pos < 33.5e6]
sig[, id := sprintf("%d:%d", chr, pos)]
# one window per locus: non-HLA signals within 2 x FLANK merged (e.g. the two CFH signals); one HLA window
setorder(sig, chr, pos)
nh <- sig[hla == FALSE]
nh[, win := cumsum(c(TRUE, diff(pos) > 2 * FLANK | diff(chr) != 0))]
win <- rbind(nh[, .(locus = paste(unique(sub(" \\(2\\)", "", locus)), collapse = "+"), chr = chr[1],
                    start = max(1L, min(pos) - FLANK), end = max(pos) + FLANK, hla = FALSE), by = win][, win := NULL],
             sig[hla == TRUE, .(locus = "HLA", chr = 6L, start = min(pos) - FLANK, end = max(pos) + FLANK, hla = TRUE)])
win[, wname := gsub("[/+ ]", "_", locus)]
fwrite(win, file.path(atl_dir, "windows_hg19.tsv"), sep = "\t")
message(nrow(sig), " signals in ", nrow(win), " windows")

## ---- Part B: one pass per GWAS -> per-window caches ---------------------------
for (a in names(gwas_files)) {
  outs <- file.path(win_dir, sprintf("%s_GWAS_%s.tsv", win$wname, a))
  if (all(file.exists(outs))) { message("cached: GWAS ", a); next }
  f <- gwas_files[[a]]; if (!file.exists(f)) stop("missing: ", f)
  hdr <- names(fread(f, nrows = 2)); ci <- match("CHR", hdr); pi <- match("BP_hg19", hdr)
  cc <- paste(sprintf("($%d==%d && $%d>=%d && $%d<=%d)", ci, win$chr, pi, win$start, pi, win$end), collapse = " || ")
  message("Scanning GWAS ", a, " (one pass) ...")
  x <- fread(cmd = sprintf("awk 'NR==1 || %s' %s", cc, shQuote(path.expand(f))))
  for (i in seq_len(nrow(win)))
    fwrite(x[CHR == win$chr[i] & BP_hg19 >= win$start[i] & BP_hg19 <= win$end[i]], outs[i], sep = "\t")
  message("  ", a, ": ", nrow(x), " rows in windows")
}
rd <- function(w, a) { g <- fread(file.path(win_dir, sprintf("%s_GWAS_%s.tsv", w, a)))
  g[, `:=`(A1 = toupper(A1), A2 = toupper(A2), id = sprintf("%d:%d", CHR, BP_hg19))]
  g[is.finite(BETA) & SE > 0] }

## ---- Part C: Aim 1 at the published lead of each signal -----------------------
align <- function(g, s) {         # effect per published risk allele
  r <- g[id == s$id]
  if (nrow(r) > 1) r <- r[toupper(A1) == s$risk | toupper(A2) == s$risk][1]
  if (!nrow(r)) return(data.table(b = NA_real_, se = NA_real_, p = NA_real_, other = NA_character_))
  sgn <- if (r$A1 == s$risk) 1 else if (r$A2 == s$risk) -1 else NA
  data.table(b = sgn * r$BETA, se = r$SE, p = r$P, other = if (r$A1 == s$risk) r$A2 else r$A1) }
a1 <- rbindlist(lapply(seq_len(nrow(sig)), function(i) {
  s <- sig[i]; w <- win[chr == s$chr & start <= s$pos & end >= s$pos, wname][1]
  e <- align(rd(w, "EUR"), s); k <- align(rd(w, "EAS"), s)
  data.table(s[, .(locus, rsid, chr, pos, id, risk, f_eur, f_eas, hla)], other = fcoalesce(e$other, k$other),
             b_eur = e$b, se_eur = e$se, p_eur = e$p, b_eas = k$b, se_eas = k$se, p_eas = k$p) }))
a1[, `:=`(or_eur = exp(b_eur), or_eas = exp(b_eas),
          q_p = pchisq((b_eur - b_eas)^2 / (se_eur^2 + se_eas^2), 1, lower.tail = FALSE))]
# contribution on the log-odds scale: V = 2p(1-p) b^2 ; exact midpoint (Shapley) split of V_EAS - V_EUR
a1[, `:=`(h_eur = 2 * f_eur * (1 - f_eur), h_eas = 2 * f_eas * (1 - f_eas))]
a1[, `:=`(V_eur = h_eur * b_eur^2, V_eas = h_eas * b_eas^2)]
a1[, `:=`(dV = V_eas - V_eur,
          d_freq   = (h_eas - h_eur) * (b_eas^2 + b_eur^2) / 2,
          d_effect = (b_eas^2 - b_eur^2) * (h_eas + h_eur) / 2)]
a1[, q_fdr := p.adjust(q_p, "BH")]
a1[, class := fifelse(is.na(dV), "missing",
              fifelse(q_fdr < 0.05, "effect differs (FDR < 0.05)",
              fifelse(abs(d_freq) >= abs(d_effect), "frequency-dominant", "effect-dominant, not significant")))]
fwrite(a1, file.path(atl_dir, "Aim1_ancestry_decomposition.tsv"), sep = "\t")
print(a1[, .(locus, rsid, OR_EUR = round(or_eur, 2), OR_EAS = round(or_eas, 2), f_eur, f_eas,
             Q_P = signif(q_p, 2), d_freq = signif(d_freq, 2), d_effect = signif(d_effect, 2), class)])
chk <- a1[!is.na(dV), sum(dV)]; message(sprintf("sum dV = %.4f ; sum(freq + effect) = %.4f (must match)", chk, a1[!is.na(dV), sum(d_freq + d_effect)]))

## ---- Part D: ancestry-specific leads (+/-250 kb of each published lead) --------
a1d <- rbindlist(lapply(seq_len(nrow(sig)), function(i) {
  s <- sig[i]; w <- win[chr == s$chr & start <= s$pos & end >= s$pos, wname][1]
  g <- list(EUR = rd(w, "EUR"), EAS = rd(w, "EAS"))
  rbindlist(lapply(names(g), function(a) {
    o <- setdiff(names(g), a); x <- g[[a]][abs(BP_hg19 - s$pos) <= 250000][order(P)][1]
    if (!nrow(x)) return(NULL)
    y <- g[[o]][id == x$id & ((A1 == x$A1 & A2 == x$A2) | (A1 == x$A2 & A2 == x$A1))][1]
    b_o <- if (nrow(y) && !is.na(y$BETA)) (if (y$A1 == x$A1) y$BETA else -y$BETA) else NA_real_
    data.table(locus = s$locus, lead_in = a, lead = x$id, A1 = x$A1, b = x$BETA, p = x$P,
               b_other = b_o, se_other = if (nrow(y)) y$SE else NA_real_, p_other = if (nrow(y)) y$P else NA_real_,
               q_p = if (!is.na(b_o)) pchisq((x$BETA - b_o)^2 / (x$SE^2 + y$SE^2), 1, lower.tail = FALSE) else NA_real_,
               dist_to_published = x$BP_hg19 - s$pos) })) }))
fwrite(a1d, file.path(atl_dir, "Aim1_ancestry_specific_leads.tsv"), sep = "\t")

## ---- Figure: frequency vs effect components (non-HLA) --------------------------
pd <- a1[hla == FALSE & !is.na(dV)]
COL <- c(`effect differs (FDR < 0.05)` = "#eb6834", `frequency-dominant` = "#2a78d6", `effect-dominant, not significant` = "#9e9e9e")
f1 <- ggplot(pd, aes(d_freq * 1e3, d_effect * 1e3, colour = class)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey60") + geom_vline(xintercept = 0, linewidth = .3, colour = "grey60") +
  geom_abline(slope = c(-1, 1), intercept = 0, linetype = 3, linewidth = .3, colour = "grey70") +
  geom_point(size = 2.2) + ggrepel::geom_text_repel(aes(label = locus), size = 2.2, colour = "grey20", fontface = "italic", max.overlaps = 30, segment.size = .2) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(x = "Frequency component of EAS - EUR contribution (x10^-3)", y = "Effect-size component (x10^-3)") +
  theme_classic(base_size = 8) + theme(legend.position = "bottom")
ggsave(file.path(atl_dir, "Aim1_decomposition.pdf"), f1, width = 120, height = 110, units = "mm")

## ---- Part E: windows in GRCh38 for ImmuNexUT / Japan Omics Browser -------------
if (!requireNamespace("rtracklayer", quietly = TRUE)) { if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("rtracklayer") }
ch_f <- file.path(data_dir, "hg19ToHg38.over.chain")
if (!file.exists(ch_f)) { download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", paste0(ch_f, ".gz"), mode = "wb")
  R.utils::gunzip(paste0(ch_f, ".gz"), remove = TRUE) }
ch <- rtracklayer::import.chain(ch_f)
gr <- GenomicRanges::GRanges(paste0("chr", win$chr), IRanges::IRanges(win$start, win$end), wname = win$wname)
l38 <- unlist(rtracklayer::liftOver(gr, ch))
w38 <- as.data.table(l38)[, .(chr38 = as.character(seqnames[1]), start38 = min(start), end38 = max(end)), by = wname]
w38 <- merge(win, w38, by = "wname", sort = FALSE)
w38[, `:=`(start38 = pmax(1L, start38 - 50000L), end38 = end38 + 50000L)]   # 50 kb margin
fwrite(w38[, .(name = wname, chr = chr38, start = start38, end = end38)], file.path(atl_dir, "windows_hg38.tsv"), sep = "\t")
message("Step 10 done: ", atl_dir)
