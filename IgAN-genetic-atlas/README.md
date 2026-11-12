# IgA nephropathy genetic atlas

Analysis code for:

> Chen H-L, *et al.* **Ancestry-resolved mapping of IgA nephropathy risk loci to immune cells and circulating proteins: APRIL, CD89 and factor H.** (submitted)

Version 1.2.0 is the version for the Genes & Immunity submission (see *Changes in v1.2.0* below).

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
| 25 | `analysis/25_matched_testing.R` | Matched-testing comparison of immune-cell and kidney colocalization |
| 26 | `analysis/26_gwas_finemap.R` | Fine-mapping of the GWAS signals |
| 27 | `analysis/27_truncation_sensitivity.R` | Truncation benchmark: complete ImmuNexUT and eQTL Catalogue data re-analyzed after keeping only variants with QTL *P* < 0.05 (Supplementary Table S12). Run after Step 11b Part 1 and the function definitions of Step 15 |
| 28 | `analysis/28_additional_analyses.py` | From Supplementary Data S1: exact prior sensitivity (p12 10⁻⁶ to 10⁻⁴), immune vs kidney counts by prior, exact matched analyses with and without the *CFH* window, parametric-bootstrap 95% CIs for the ancestry decomposition (20 000 replicates), Steiger directionality |
| 29 | `analysis/29_supplementary_data_S1.py` | Builds Supplementary Data S1 (Tables S1–S10 and S12) |
| Figs | `figures/20_fig1_design.R` … `figures/22_fig6_complement.R` | Figures 1–6 (`21_fig2_ancestry.R` now reads the bootstrap CIs from Step 28) |
| Fig S2 | `figures/29_figS2_fcar_transcript.py` | *FCAR* transcript structure and the rs1865097 junction (Ensembl ENST00000355524, GRCh38) |

Notes
- AlphaGenome splicing scores are unsigned magnitudes; activity (`*_ACTIVE`) scorers describe predicted activity, not allelic effect, and are not re-signed. The per-variant summary printed by Step 14 is superseded by the summaries in the manuscript and Supplemental Table S7.
- Positions are GRCh37 unless stated; GRCh38 resources are lifted with UCSC liftOver (`hg38ToHg19.over.chain`).

## Results files

| File | Content |
|---|---|
| `results/Supplementary_Tables_S1-S8_complete.xlsx` | Complete tables from v1.1.3 |
| `results/Supplementary_Data_S1_v1.1.3.xlsx` | Supplementary Data S1 as first submitted (input to Steps 28–29) |
| `results/Supplementary_Data_S1.xlsx` | Supplementary Data S1 of the Genes & Immunity submission (Tables S1–S10 and S12) |
| `results/additional_analyses.xlsx` | All Step 28 output tables |
| `results/Aim1_decomposition_bootstrapCI.tsv` | Bootstrap 95% CIs for the ancestry decomposition (Figure 2) |
| `results/truncation/` | Step 27 per-test output (`truncation_immunexut.tsv`, `truncation_eqtlcatalogue.tsv`) |

## Changes in v1.2.0

- Added Steps 27–29 and the Supplementary Figure S2 script.
- Figure 2 shows parametric-bootstrap 95% CIs.
- Prior sensitivity extended to p12 = 10⁻⁴; matched analyses repeated without the whole *CFH* window.
- Truncation benchmark for the Japan COVID-19 Task Force e/pQTL data (truncated at *P* < 0.05).
- No change to the v1.1.3 primary results.

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
Python 3.13 (alphagenome, remotezip, pyliftover, pandas, numpy, scipy, openpyxl, matplotlib).

## License

MIT (see `LICENSE`).
