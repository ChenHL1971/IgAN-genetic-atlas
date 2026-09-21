# =============================================================================
# Run all Kidney International figures in one go (Figures 1-6 + Supplementary Figure S1).
# Usage in a fresh R session:  source("~/Downloads/figure_scripts_KI/run_KI_figures.R")
# Output: ~/Desktop/IgAN_Genetics_Paper2/results/figures/  (PDF + 600-dpi TIFF for each figure)
# =============================================================================
scripts_dir <- "~/Downloads/figure_scripts_KI"          # folder created when you double-click the zip
if (!dir.exists(path.expand(scripts_dir))) stop("Folder not found: ", scripts_dir, " - change scripts_dir on line 6")

## ---- shared settings and helper functions (from Step 11b, Parts A-B; uses cached files, no downloads) ----
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

## ---- run each figure; one failure does not stop the others ----
figs <- c("20_fig1_design.R", "21_fig2_ancestry.R", "17_fig3_atlas_map.R", "18_fig4_april_chain.R",
          "19_fig5_fcar.R", "22_fig6_complement.R", "25_matched_testing.R")
status <- sapply(figs, function(f) {
  message("\n==== ", f, " ====")
  tryCatch({ source(file.path(path.expand(scripts_dir), f), local = new.env(parent = globalenv())); "OK" },
           error = function(e) paste("ERROR:", conditionMessage(e)))
})
cat("\n==== Summary ====\n"); print(data.frame(script = figs, status = status, row.names = NULL))
cat("\nFigures are in:", path.expand(file.path(proj_dir, "results", "figures")), "\n")
if (capabilities("aqua")) system(paste("open", shQuote(path.expand(file.path(proj_dir, "results", "figures")))))
