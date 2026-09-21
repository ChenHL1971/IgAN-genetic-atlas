## ---- Fig. 3: cell-type and tissue map of IgAN risk loci (paste after Step 11b Part 1; needs atl_dir, proj_dir) ----
for (p in c("data.table", "ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table); library(ggplot2); library(patchwork)
FAM <- if (capabilities("aqua")) "Arial" else "sans"   # Kidney International: Arial; lowercase panel labels
fig_dir <- file.path(proj_dir, "results", "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
SEQ <- c("#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")

## 1. read the four result tables -> one long table: window, ancestry, group, column, gene, PP.H4
imx <- fread(file.path(atl_dir, "Atlas_immune_coloc_all.tsv.gz"))
kid <- fread(file.path(atl_dir, "Atlas_kidney_coloc_all.tsv.gz"))
ecr <- fread(file.path(atl_dir, "Atlas_eqtlcatalogue_coloc.tsv.gz"))
job <- fread(file.path(atl_dir, "Atlas_JOB_pqtl_eqtl.tsv"))

imx_col <- c(CL_Mono = "Monocytes", Int_Mono = "Monocytes", NC_Mono = "Monocytes", CD16p_Mono = "Monocytes",
             mDC = "mDC", pDC = "pDC", Neu = "Neutrophils", LDG = "Neutrophils",
             Naive_B = "B cells", USM_B = "B cells", SM_B = "B cells", DN_B = "B cells", Plasmablast = "Plasmablasts",
             Naive_CD4 = "CD4 T", Mem_CD4 = "CD4 T", Th1 = "CD4 T", Th2 = "CD4 T", Th17 = "CD4 T", Tfh = "CD4 T",
             Fr_I_nTreg = "Treg", Fr_II_eTreg = "Treg", Fr_III_T = "Treg",
             Naive_CD8 = "CD8 T", CM_CD8 = "CD8 T", EM_CD8 = "CD8 T", TEMRA_CD8 = "CD8 T", Mem_CD8 = "CD8 T", NK = "NK")
L1 <- imx[, .(window, ancestry, group = "Japanese immune cells\n(ImmuNexUT)", column = imx_col[cell], gene, PP.H4)]

ec_col <- function(tissue, quant) fcase(
  quant == "leafcutter" & tissue == "liver", "Liver",
  quant == "leafcutter" & tissue == "kidney cortex", "Cortex\n(GTEx)",
  quant == "leafcutter", "Splicing\n(immune, blood)",
  tissue %in% c("monocyte", "CD16+ monocyte"), "Monocytes",
  tissue == "neutrophil", "Neutrophils",
  tissue %in% c("B cell", "memory B cell"), "B cells",
  tissue == "CD4+ T cell", "CD4 T",
  tissue == "NK cell", "NK",
  tissue %in% c("dendritic cell", "plasmacytoid dendritic cell"), "DC",
  tissue == "blood", "Whole blood",
  tissue == "liver", "Liver",
  tissue == "kidney cortex", "Cortex\n(GTEx)", default = NA_character_)
L2 <- ecr[, .(window, ancestry, column = ec_col(tissue, quant), gene, PP.H4)]
L2[, group := fifelse(column == "Cortex\n(GTEx)", "Kidney", "European cells and tissues\n(eQTL Catalogue)")]
L3 <- job[, .(window, ancestry, group = "Japanese\nblood\n(JCTF)", column = fifelse(layer == "pQTL", "Plasma\nprotein", "Whole-blood\nmRNA"), gene, PP.H4)]
L4 <- kid[, .(window, ancestry, group = "Kidney", column = fifelse(tissue == "Glom", "Glomerulus", "Tubulo-\ninterstitium"), gene, PP.H4)]
A <- rbindlist(list(L1, L2, L3, L4), use.names = TRUE)[!is.na(column) & !is.na(PP.H4)]

## 2. best gene per window x ancestry x column
A[grepl("^ENSG", gene), gene := "(lncRNA)"]
B <- A[order(-PP.H4)][, .SD[1], by = .(window, ancestry, group, column)]

## 3. axis order
col_order <- c("Monocytes", "mDC", "pDC", "Neutrophils", "B cells", "Plasmablasts", "CD4 T", "Treg", "CD8 T", "NK",
               "DC", "Whole blood", "Splicing\n(immune, blood)", "Liver",
               "Whole-blood\nmRNA", "Plasma\nprotein",
               "Glomerulus", "Tubulo-\ninterstitium", "Cortex\n(GTEx)")
grp_order <- c("Japanese immune cells\n(ImmuNexUT)", "European cells and tissues\n(eQTL Catalogue)", "Japanese\nblood\n(JCTF)", "Kidney")
row_order <- c("TNFSF12_13", "CARD9", "ITGAM_ITGAX", "FCAR", "LIF_OSM", "TNFSF8_15", "PF4V1", "IRF8",
               "LYN", "FCRL3", "TNFSF4_18", "CD28", "REL", "ZMIZ1", "IRF4_DUSP22", "OVOL1_RELA", "REEP3",
               "TNFRSF13B", "CFH", "IGH", "ANXA3", "DEFA1_4", "ETS1", "LY86")
row_lab <- c(TNFSF12_13 = "TNFSF12/13", ITGAM_ITGAX = "ITGAM/ITGAX", LIF_OSM = "LIF/OSM", TNFSF8_15 = "TNFSF8/15",
             TNFSF4_18 = "TNFSF4/18", IRF4_DUSP22 = "IRF4/DUSP22", OVOL1_RELA = "OVOL1/RELA", DEFA1_4 = "DEFA1/4")
# ImmuNexUT and eQTL Catalogue share some column names, so the grid is built per group
eu_cols <- c("Monocytes", "Neutrophils", "B cells", "CD4 T", "NK", "DC", "Whole blood", "Splicing\n(immune, blood)", "Liver")
full <- rbind(CJ(window = row_order, ancestry = c("EAS", "EUR"), column = col_order[1:10])[, group := grp_order[1]],
              CJ(window = row_order, ancestry = c("EAS", "EUR"), column = eu_cols)[, group := grp_order[2]],
              CJ(window = row_order, ancestry = c("EAS", "EUR"), column = col_order[15:16])[, group := grp_order[3]],
              CJ(window = row_order, ancestry = c("EAS", "EUR"), column = col_order[17:19])[, group := grp_order[4]])
D <- merge(full, B, by = c("window", "ancestry", "group", "column"), all.x = TRUE)
D[, `:=`(window = factor(window, rev(row_order)), group = factor(group, grp_order),
         column = factor(column, unique(c(col_order[1:10], eu_cols, col_order[15:19]))),
         anc_lab = factor(fifelse(ancestry == "EAS", "East Asian GWAS", "European GWAS"), c("East Asian GWAS", "European GWAS")))]
D[, txtcol := fifelse(!is.na(PP.H4) & PP.H4 >= 0.9, "white", "grey10")]
## genes behind the strong cells, listed to the right of each row (most frequent first; lncRNAs omitted)
G <- D[PP.H4 >= 0.8 & gene != "(lncRNA)", .N, by = .(window, anc_lab, gene)][order(-N)]
G <- G[, .(genes = paste(head(gene, 3), collapse = ", ")), by = .(window, anc_lab)]
G <- merge(unique(D[, .(window, anc_lab)]), G, by = c("window", "anc_lab"), all.x = TRUE)[is.na(genes), genes := ""]

## 4. plot
p3 <- ggplot(D, aes(column, window)) +
  geom_tile(aes(fill = PP.H4), colour = "white", linewidth = 0.4) +
  geom_point(data = D[PP.H4 >= 0.8], aes(colour = txtcol), size = 0.6) +
  scale_colour_identity() +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 1), breaks = c(0, 0.5, 0.8, 1), na.value = "#e6e6e6", name = "PP.H4") +
  scale_y_discrete(labels = function(x) ifelse(x %in% names(row_lab), row_lab[x], x)) +
  facet_grid(anc_lab ~ group, scales = "free_x", space = "free_x", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7, base_family = FAM) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(face = "italic", size = 6), strip.placement = "outside",
        strip.text.x = element_text(face = "bold", size = 6.5), strip.text.y.left = element_text(face = "bold", size = 7, angle = 90),
        panel.spacing.x = unit(1.2, "mm"), panel.spacing.y = unit(3, "mm"),
        legend.position = "bottom", legend.key.height = unit(2.5, "mm"), legend.key.width = unit(10, "mm"))
pg <- ggplot(G, aes(0, window, label = genes)) + geom_text(hjust = 0, size = 1.9, fontface = "italic", colour = "grey15") +
  facet_grid(anc_lab ~ "Genes\n(PP.H4 \u2265 0.8)") + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_size = 7) + theme(strip.text.x = element_text(face = "bold", size = 6.5, vjust = 0), strip.text.y = element_blank(),
                                     panel.spacing.y = unit(3, "mm"))
p3 <- p3 + pg + plot_layout(widths = c(5, 1.3))
if (capabilities("aqua")) { quartz(type = "pdf", file = file.path(fig_dir, "Fig3_atlas_map.pdf"), width = 180/25.4, height = 235/25.4); print(p3); dev.off()
} else ggsave(file.path(fig_dir, "Fig3_atlas_map.pdf"), p3, width = 180, height = 235, units = "mm", bg = "white", device = cairo_pdf)
ggsave(file.path(fig_dir, "Fig3_atlas_map.tiff"), p3, width = 180, height = 235, units = "mm", dpi = 600, compression = "lzw", bg = "white")
fwrite(D[, .(window, ancestry, group, column = gsub("\n", " ", column), gene, PP.H4)], file.path(fig_dir, "Fig3_source_data.tsv"), sep = "\t")
cat("Fig 3 saved to", fig_dir, "\n")
cat("Loci with any PP.H4 >= 0.8, by group:\n"); print(D[PP.H4 >= 0.8, .(loci = uniqueN(window)), by = group])
