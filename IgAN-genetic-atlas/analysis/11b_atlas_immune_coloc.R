# =============================================================================
# Paper 2 (atlas) - Step 11b: IgAN GWAS (EUR, EAS) x ImmuNexUT (28 cell types)
#   for every non-HLA atlas window from Step 10.
#   Part A  1000G phase 3 genotypes per window (remote tabix; no whole-chromosome download)
#           -> compact .rds: SNV positions/alleles + 0/1/2 dosages for EUR, EAS (incl. JPT)
#   Part B  split ImmuNexUT atlas extraction (Step 11a) by window
#   Part C  coloc.abf per window x ancestry x cell x gene (>= 50 shared SNPs)
#   Part D  summary tables (+ LIF/OSM, TNFSF12/13, FCAR, LYN first look)
# Needs Step 10 outputs (results/atlas/windows_hg19.tsv, windows/<w>_GWAS_<anc>.tsv)
# =============================================================================
proj_dir <- "~/Desktop/IgAN_Genetics_Paper2"
atl_dir  <- file.path(proj_dir, "results", "atlas")
win_dir  <- file.path(atl_dir, "windows")
imx_in   <- file.path(proj_dir, "data", "ImmuNexUT", "atlas_hg38")
kg_dir   <- file.path(atl_dir, "kg"); dir.create(kg_dir, showWarnings = FALSE)
co_dir   <- file.path(atl_dir, "coloc_immune"); dir.create(co_dir, showWarnings = FALSE)
panel_f  <- file.path(proj_dir, "data", "1000G", "integrated_call_samples_v3.20130502.ALL.panel")
chain_f  <- file.path(proj_dir, "data", "hg38ToHg19.over.chain")
N_imx <- 416; ncase <- c(EUR = 5556, EAS = 4590); nctrl <- c(EUR = 21178, EAS = 7573)
KG_URL <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

for (p in c("data.table", "coloc", "BiocManager")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in c("Rsamtools", "rtracklayer")) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
library(data.table); library(coloc)

win <- fread(file.path(atl_dir, "windows_hg19.tsv"))[hla == FALSE]
panel <- fread(panel_f, header = FALSE, fill = TRUE, select = 1:4, col.names = c("sample", "pop", "super_pop", "gender"))[sample != "sample"]
keep_s <- panel[super_pop %in% c("EUR", "EAS"), sample]

## ---- Part A: 1000G per window (remote tabix) ------------------------------------
kg_window <- function(w) {
  f <- file.path(kg_dir, paste0(w, ".rds")); if (file.exists(f)) return(readRDS(f))
  r <- win[wname == w]
  tf <- Rsamtools::TabixFile(sprintf(KG_URL, r$chr))
  hdr <- Rsamtools::headerTabix(tf)$header; cols <- strsplit(tail(hdr, 1), "\t")[[1]]
  gr <- GenomicRanges::GRanges(as.character(r$chr), IRanges::IRanges(r$start, r$end))
  message(sprintf("  1000G %s: fetching chr%s:%d-%d ...", w, r$chr, r$start, r$end))
  ln <- Rsamtools::scanTabix(tf, param = gr)[[1]]
  x <- fread(text = c(paste(cols, collapse = "\t"), ln), sep = "\t", header = TRUE,
             select = c("#CHROM", "POS", "REF", "ALT", intersect(keep_s, cols)), colClasses = list(character = c("REF", "ALT")))
  setnames(x, "#CHROM", "CHROM")
  x <- x[nchar(REF) == 1 & nchar(ALT) == 1][!duplicated(POS) & !duplicated(POS, fromLast = TRUE)]
  S <- intersect(keep_s, names(x)); g <- as.matrix(x[, ..S])
  D <- (substr(g, 1, 1) == "1") + (substr(g, 3, 3) == "1"); storage.mode(D) <- "integer"; colnames(D) <- S
  out <- list(meta = data.table(SNP = paste(x$CHROM, x$POS, sep = ":"), POS19 = x$POS, kREF = x$REF, kALT = x$ALT), D = D)
  saveRDS(out, f); message(sprintf("  1000G %s: %d SNVs x %d samples saved", w, nrow(D), ncol(D))); out
}
kg_freq <- function(k) {
  s <- colnames(k$D); jpt <- s %in% panel[pop == "JPT", sample]; eas <- s %in% panel[super_pop == "EAS", sample]
  k$meta[, `:=`(AF_JPT = rowMeans(k$D[, jpt, drop = FALSE]) / 2, AF_EAS = rowMeans(k$D[, eas, drop = FALSE]) / 2,
                AF_EUR = rowMeans(k$D[, !eas, drop = FALSE]) / 2)][]
}

## ---- Part B: split ImmuNexUT by window --------------------------------------------
imx_split <- function() {
  done <- file.path(co_dir, "_split_done"); if (file.exists(done)) return(invisible(message("cached: ImmuNexUT split")))
  fs <- list.files(imx_in, pattern = "\\.tsv$", full.names = TRUE); if (!length(fs)) stop("no ImmuNexUT files in ", imx_in, " (run Step 11a)")
  message(sprintf("Reading %d cell-type files ...", length(fs)))
  x <- rbindlist(lapply(fs, function(f) fread(f)[, cell := sub("\\.tsv$", "", basename(f))]))
  setnames(x, make.names(names(x)))
  for (w in unique(x$window)) fwrite(x[window == w], file.path(co_dir, paste0("imx_", w, "_hg38.tsv")), sep = "\t")
  print(x[, .N, by = window]); file.create(done)
}
ch <- NULL
lift19 <- function(chr, pos) {
  if (is.null(ch)) ch <<- rtracklayer::import.chain(path.expand(chain_f))
  l <- rtracklayer::liftOver(GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos)), ch); n <- lengths(l)
  out <- rep(NA_integer_, length(pos)); out[n == 1] <- as.integer(BiocGenerics::start(unlist(l[n == 1]))); out
}
get_gwas <- function(w, a) {
  x <- fread(file.path(win_dir, sprintf("%s_GWAS_%s.tsv", w, a)))
  x <- x[!is.na(BETA) & SE > 0 & !is.na(P)][, `:=`(A1 = toupper(A1), A2 = toupper(A2), SNP = paste(CHR, BP_hg19, sep = ":"))]
  x <- x[!((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A") | (A1 == "C" & A2 == "G") | (A1 == "G" & A2 == "C"))]
  unique(x, by = "SNP")[, varbeta := SE^2]
}

## ---- Part C: coloc per window ----------------------------------------------------------
imx_window <- function(w, min_snps = 50) {
  of <- file.path(co_dir, paste0("coloc_", w, ".tsv")); if (file.exists(of)) return(fread(of))
  ef <- file.path(co_dir, paste0("imx_", w, "_hg38.tsv")); if (!file.exists(ef)) { message("no ImmuNexUT rows: ", w); return(NULL) }
  e <- fread(ef)
  e <- e[!is.na(slope.ALT.) & nominal_P_value > 0 & nominal_P_value < 1 & nchar(REF) == 1 & nchar(ALT) == 1]
  e[, z := qnorm(nominal_P_value / 2, lower.tail = FALSE) * sign(slope.ALT.)]; e <- e[abs(z) > 0][, se := abs(slope.ALT. / z)]
  up <- unique(e[, .(Variant_CHR, Variant_position_start)])[, pos19 := lift19(Variant_CHR, Variant_position_start)]
  e <- merge(e, up[!is.na(pos19)], by = c("Variant_CHR", "Variant_position_start"))
  e[, SNP := paste(sub("^chr", "", Variant_CHR), pos19, sep = ":")]
  e <- merge(e, kg_freq(kg_window(w)), by = "SNP")
  e[, kdir := fifelse(kREF == REF & kALT == ALT, 1, fifelse(kREF == ALT & kALT == REF, -1, NA_real_))]
  e <- e[!is.na(kdir)][, af_alt := fifelse(kdir == 1, AF_JPT, 1 - AF_JPT)][, MAF := pmin(af_alt, 1 - af_alt)][MAF >= 0.01]
  e <- unique(e, by = c("cell", "Gene_id", "SNP"))
  res <- list()
  for (a in c("EAS", "EUR")) {
    g <- get_gwas(w, a)
    m <- merge(e, g[, .(SNP, gA1 = A1, gA2 = A2, BETA, varbeta, P)], by = "SNP")
    m[, gdir := fifelse(gA1 == ALT & gA2 == REF, 1, fifelse(gA1 == REF & gA2 == ALT, -1, NA_real_))]
    m <- m[!is.na(gdir)][, gwas_beta := BETA * gdir]
    grp <- m[, .N, by = .(cell, Gene_id, Gene_name)][N >= min_snps]
    for (i in seq_len(nrow(grp))) {
      d <- m[cell == grp$cell[i] & Gene_id == grp$Gene_id[i]]
      invisible(capture.output(r <- suppressMessages(suppressWarnings(coloc.abf(
        list(beta = d$gwas_beta, varbeta = d$varbeta, snp = d$SNP, type = "cc", s = ncase[[a]] / (ncase[[a]] + nctrl[[a]])),
        list(beta = d$slope.ALT., varbeta = d$se^2, snp = d$SNP, type = "quant", N = N_imx, MAF = d$MAF))))))
      sm <- r$summary; top <- r$results[which.max(r$results$SNP.PP.H4), ]; ti <- match(top$snp, d$SNP)
      res[[length(res) + 1]] <- data.table(window = w, ancestry = a, cell = grp$cell[i], gene = grp$Gene_name[i], gene_id = grp$Gene_id[i],
        nsnps = nrow(d), min_gwas_p = min(d$P), min_eqtl_p = min(d$nominal_P_value),
        PP.H3 = sm[["PP.H3.abf"]], PP.H4 = sm[["PP.H4.abf"]], top_snp_H4 = top$snp, top_SNP.PP.H4 = top$SNP.PP.H4,
        dir_at_top = sign(d$slope.ALT.[ti]) * sign(d$gwas_beta[ti]))
    }
    message(sprintf("  %s %s: %d cell x gene pairs", w, a, nrow(grp)))
  }
  res <- rbindlist(res); if (nrow(res)) fwrite(res, of, sep = "\t"); res
}

## ---- run --------------------------------------------------------------------------------
imx_split()
t0 <- Sys.time()
allr <- rbindlist(lapply(win$wname, function(w) tryCatch(imx_window(w), error = function(e) { message("ERROR ", w, ": ", conditionMessage(e)); NULL })), fill = TRUE)
message(sprintf("coloc done: %d tests in %.0f min", nrow(allr), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (!nrow(allr)) stop("no colocalization tests ran - check that Step 11a finished and allele matching worked")
fwrite(allr, file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"), sep = "\t")

## ---- Part D: summaries ---------------------------------------------------------------------
lineage <- c(CL_Mono = "Myeloid", Int_Mono = "Myeloid", CD16p_Mono = "Myeloid", NC_Mono = "Myeloid", mDC = "Myeloid", pDC = "Myeloid",
             Neu = "Myeloid", LDG = "Myeloid", Naive_B = "B", USM_B = "B", SM_B = "B", DN_B = "B", Plasmablast = "B", NK = "NK")
allr[, lineage := fifelse(cell %in% names(lineage), lineage[cell], "T")]
best <- allr[order(-PP.H4)][, .SD[1], by = .(window, ancestry)]
cat("\n[1] Best immune colocalization per window (PP.H4 >= 0.5 marked):\n")
print(dcast(best, window ~ ancestry, value.var = c("PP.H4", "gene", "cell"))[order(-pmax(PP.H4_EAS, PP.H4_EUR, na.rm = TRUE))], nrows = 50)
cat("\n[2] Windows x lineage, max PP.H4 (EAS | EUR):\n")
lw <- allr[, .(PP = round(max(PP.H4), 2)), by = .(window, ancestry, lineage)]
print(dcast(lw, window ~ ancestry + lineage, value.var = "PP"), nrows = 50)
for (w in c("LIF_OSM", "TNFSF12_13", "FCAR", "LYN")) if (w %in% allr$window) {
  cat("\n[3]", w, ": pairs with PP.H4 >= 0.5\n")
  print(allr[window == w & PP.H4 >= 0.5, .(ancestry, cell, gene, nsnps, min_eqtl_p = signif(min_eqtl_p, 2), PP.H3 = round(PP.H3, 2),
        PP.H4 = round(PP.H4, 2), top_snp_H4, top_pp = round(top_SNP.PP.H4, 2), dir_at_top)][order(-PP.H4)], nrows = 40)
}
message("Step 11b done: ", atl_dir)
