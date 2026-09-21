# IgA nephropathy genetic atlas

Analysis code for:

> Chen H-L, *et al.* **An ancestry-resolved cell-type and molecular map of genetic risk for IgA nephropathy.** (submitted)

All analyses use publicly available summary statistics; no individual-level data are included or required.

## Pipeline

Run in order. R scripts assume an R session in which `proj_dir` points to the project folder
(default `~/Desktop/IgAN_Genetics_Paper2`); later steps reuse objects defined in Step 11b, Part 1.

| Step | Script | Purpose |
|---|---|---|
| 10 | `analysis/10_atlas_aim1.R` | 30 IgAN signals (Kiryluk 2023); ±500 kb windows; per-window GWAS caches; ancestry decomposition (2p(1−p)β², frequency vs effect split; z test, BH-FDR) |
| 11a | `analysis/11a_immunexut_stream_atlas.py` | Stream ImmuNexUT nominal eQTL (E-GEAD-420) and keep variants in atlas windows (GRCh38) |
| 11b | `analysis/11b_atlas_immune_coloc.R` | 1000 Genomes genotypes by remote tabix; liftOver; coloc.abf with 28 immune-cell types (Part 1 defines shared helpers) |
| 12 | `analysis/12_atlas_susie.R` | SuSiE-based colocalization (coloc.susie, L = 5) for prioritized pairs |
| 13 | `analysis/13_atlas_kidney_coloc.R` | coloc.abf with NephQTL2 glomerulus and tubulointerstitium |
| 14 | `analysis/14_alphagenome_scoring.py` | AlphaGenome variant scoring (requires an AlphaGenome API key in `ALPHAGENOME_API_KEY`) |
| 15 | `analysis/15_eqtlcat_replication_sqtl.R` | eQTL Catalogue expression and LeafCutter splicing QTL (25 datasets, remote tabix) |
| 16 | `analysis/16_job_pqtl.R` | Japan COVID-19 Task Force whole-blood eQTL and Olink pQTL; coloc and cis-MR (Wald ratio) |
| — | `analysis/23_supplemental_tables.R` | Supplemental Tables S1–S7 |
| Figs | `figures/20_fig1_design.R` … `figures/22_fig6_complement.R` | Figures 1–6 |

Notes
- AlphaGenome splicing scores are unsigned magnitudes; activity (`*_ACTIVE`) scorers describe predicted activity, not allelic effect, and are not re-signed. The per-variant summary printed by Step 14 is superseded by the summaries in the manuscript and Supplemental Table S7.
- Positions are GRCh37 unless stated; GRCh38 resources are lifted with UCSC liftOver (`hg38ToHg19.over.chain`).

## Data sources

| Resource | Access |
|---|---|
| IgAN GWAS (European and East Asian meta-analyses; Kiryluk *et al.* 2023) | Kiryluk laboratory, Columbia University |
| ImmuNexUT nominal cis-eQTL | NBDC / DDBJ E-GEAD-420 |
| eQTL Catalogue (OneK1K, BLUEPRINT, Quach, Schmiedel, GTEx) | https://www.ebi.ac.uk/eqtl/ |
| Japan COVID-19 Task Force e/pQTL (Wang *et al.* 2024) | NBDC hum0343.v3.qtl.v1 |
| NephQTL2 | Han *et al.* 2023 |
| 1000 Genomes phase 3 | IGSR |
| AlphaGenome | https://github.com/google-deepmind/alphagenome |

## Software

R 4.6.1 (data.table, coloc 5.2.3, susieR, rtracklayer, Rsamtools, ggplot2, patchwork, ggrepel, openxlsx, qpdf);
Python 3.13 (alphagenome, remotezip, pyliftover, pandas).

## License

MIT (see `LICENSE`).
