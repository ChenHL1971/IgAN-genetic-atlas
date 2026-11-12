"""Step 29. Builds Supplementary Data S1 (Tables S1-S10 and S12) for the Genes & Immunity submission.
Usage: python 29_supplementary_data_S1.py Supplementary_Data_S1_v1.1.3.xlsx additional_analyses.xlsx Supplementary_Data_S1.xlsx
S12 values are the summaries of results/atlas/truncation/*.tsv from Step 27."""
import sys
import openpyxl, pandas as pd
from copy import copy
wb = openpyxl.load_workbook(sys.argv[1])
r = pd.read_excel(sys.argv[2], sheet_name=None)
hf = copy(wb['S1'].cell(1,1).font)
def write_df(ws, df, row=1, col=1):
    for j,c in enumerate(df.columns): ws.cell(row, col+j, str(c)).font = copy(hf)
    for i,rr in enumerate(df.itertuples(index=False)):
        for j,v in enumerate(rr): ws.cell(row+1+i, col+j, None if (isinstance(v,float) and pd.isna(v)) else v)
# S6 CIs
s6 = wb['S6']; ci = r['T4_decomp_per_signal']; last = s6.max_column
cols = ['V_EUR 95% CI','V_EAS 95% CI','Difference 95% CI','Frequency component 95% CI','Effect component 95% CI']
for j,c in enumerate(cols):
    s6.cell(1,last+1+j,c.replace(' ','_')).font = copy(hf)
    for i,v in enumerate(ci[c]): s6.cell(2+i,last+1+j,v)
r0 = s6.max_row+2
s6.cell(r0,1,'Totals, 25 non-HLA signals (parametric bootstrap, 20 000 replicates; allele frequencies resampled from 1000 Genomes; signals and ancestries treated as independent; winner\'s curse not corrected)').font = copy(hf)
tot = r['T4b_decomp_totals']; tot = tot.rename(columns={tot.columns[0]:'Quantity'})
tot = tot[tot.Quantity != 'Frequency share of difference (%)']
extra = pd.DataFrame({'Quantity':['Frequency component (sum)','Effect-size component (sum)','Difference V_EAS - V_EUR (sum)'],
   'Estimate':[0.04444,0.04116,0.08559],'2.5%':[0.02710,-0.03550,0.02788],'97.5%':[0.08721,0.09410,0.13902]})
extra.columns = tot.columns
write_df(s6, pd.concat([tot, extra]), row=r0+1)
# S8
s8 = wb['S8']; r0 = s8.max_row+2
s8.cell(r0,1,'Exact versions (enumeration / hypergeometric): all genes, excluding immune-cell CFH gene tests, and excluding the whole CFH window').font = copy(hf)
m = r['T3_matched_exact'].rename(columns={'Unnamed: 0':'Analysis'}); write_df(s8, m, row=r0+1)
write_df(s8, r['T3b_gene_level'], row=s8.max_row+2)
# S9
ws = wb.create_sheet('S9'); ws.cell(1,1,'Sensitivity of PP.H4 to the prior probability of a shared variant (p12); published analysis used p12 = 1e-5').font = copy(hf)
write_df(ws, r['T1_prior_main'], row=3); write_df(ws, r['T2_counts_by_prior'], row=ws.max_row+3)
# S10
ws = wb.create_sheet('S10'); ws.cell(1,1,'Steiger directionality test for single-instrument Wald-ratio estimates').font = copy(hf)
write_df(ws, r['T5_steiger'], row=3)
# S11 truncation benchmark (from Step 27 output)
ws = wb.create_sheet('S12')
ws.cell(1,1,'Truncation benchmark: complete QTL data re-analyzed after keeping only variants with QTL P < 0.05 (partial emulation of the Japan COVID-19 Task Force release rule; PIP not available). Tests with < 30 retained variants were not evaluable.').font = copy(hf)
A = pd.DataFrame([dict(tests=2414, not_evaluable=700, spearman=0.757, median_abs_diff=0.015, strong_kept=195, strong_lost=0, strong_gained=21, same_top_snp=0.516)])
ws.cell(3,1,'A. ImmuNexUT (all tests with published PP.H4 >= 0.5 plus 2 000 random tests below 0.5)').font = copy(hf); write_df(ws, A, row=4)
G = pd.DataFrame([
 ("FCRL5","EAS","DN_B",0.983,0.983,2210,688,True),("TNFSF4","EAS","Fr_III_T",0.014,0.008,1231,267,True),("TNFSF4","EUR","Naive_B",0.852,0.845,1192,199,True),
 ("CFH","EAS","NC_Mono",0.939,0.940,1551,233,True),("CFH","EUR","NC_Mono",0.548,0.926,951,136,True),("CD28","EAS","Naive_CD4",0.981,0.981,1366,317,True),
 ("CD28","EUR","Naive_CD4",0.702,0.709,1256,308,True),("LYN","EAS","Plasmablast",0.993,0.993,1891,643,True),("LYN","EUR","Plasmablast",0.998,0.998,1750,629,True),
 ("CARD9","EAS","Mem_CD8",0.009,0.001,1789,78,False),("CARD9","EUR","NC_Mono",0.974,0.974,1694,474,True),("ITGAM","EAS","CL_Mono",0.935,0.935,1070,160,True),
 ("ITGAM","EUR","CL_Mono",0.932,0.932,1038,156,True),("IRF8","EAS","pDC",0.767,0.783,2733,243,True),("IRF8","EUR","pDC",0.908,0.930,2486,233,True),
 ("TNFSF13","EAS","CL_Mono",1.000,1.000,1657,467,True),("TNFSF13","EUR","EM_CD8",0.400,0.414,1569,433,True)],
 columns=["gene","GWAS","cell","PP.H4_full","PP.H4_truncated","variants_full","variants_truncated","same_top_variant"])
ws.cell(ws.max_row+2,1,'A2. Prioritized genes (best cell type)').font = copy(hf); write_df(ws, G, row=ws.max_row+1)
C = pd.DataFrame([("REL","EUR","Naive_CD4","AC016747.1",0.796,0.825,396),("TNFRSF13B","EUR","Naive_B","SREBF1",0.795,0.812,140),("LYN","EUR","Plasmablast","AC018607.1",0.794,0.820,138),
 ("CD28","EAS","Mem_CD4","CD28",0.792,0.800,141),("TNFRSF13B","EAS","Th17","ZNF624",0.791,0.858,281),("TNFRSF13B","EAS","NK","ZNF624",0.790,0.858,285),("ZMIZ1","EUR","Mem_CD4","PPIF",0.789,0.809,186),
 ("TNFRSF13B","EAS","Tfh","ZNF624",0.788,0.857,229),("TNFRSF13B","EAS","Naive_CD8","ZNF624",0.783,0.853,354),("TNFRSF13B","EAS","Th2","ZNF624",0.783,0.854,138),("TNFRSF13B","EAS","Naive_CD4","ZNF624",0.782,0.853,147),
 ("TNFRSF13B","EAS","DN_B","ZNF624",0.780,0.851,143),("TNFRSF13B","EAS","Fr_III_T","ZNF624",0.780,0.851,142),("TNFRSF13B","EAS","Fr_I_nTreg","ZNF624",0.776,0.848,236),("TNFRSF13B","EAS","Naive_B","ZNF624",0.771,0.844,148),
 ("TNFRSF13B","EAS","Mem_CD4","ZNF624",0.757,0.834,245),("TNFRSF13B","EAS","Plasmablast","ZNF624",0.745,0.808,229),("TNFRSF13B","EAS","USM_B","ZNF624",0.741,0.822,203),("LY86","EAS","Neu","RREB1",0.712,0.840,304),
 ("CFH","EUR","NC_Mono","CFH",0.548,0.926,136),("CARD9","EUR","NK","CARD9",0.035,0.967,449)],
 columns=["window","GWAS","cell","gene","PP.H4_full","PP.H4_truncated","variants_truncated"])
ws.cell(ws.max_row+2,1,'A3. Tests that crossed PP.H4 0.8 only after truncation (21)').font = copy(hf); write_df(ws, C, row=ws.max_row+1)
B = pd.DataFrame([("GTEx blood TNFSF13","EAS",1552,317,0.001,0.021,False),("GTEx blood TNFSF13","EUR",2163,388,0.158,0.171,True),("GTEx blood CFH","EAS",1233,40,0.044,0.114,False),
 ("GTEx blood CFH","EUR",1463,57,0.041,0.003,False),("GTEx liver CFHR3","EAS",1170,405,0.821,0.823,True),("GTEx liver CFHR3","EUR",1459,568,0.601,0.592,True),
 ("GTEx liver CFHR4","EAS",1170,565,0.282,0.341,True),("GTEx liver CFHR4","EUR",1459,632,0.828,0.828,True),("Quach monocyte FCAR junction","EAS",2272,191,0.969,0.987,True),
 ("Quach monocyte FCAR junction","EUR",3167,252,0.681,0.725,True),("GTEx blood FCAR junction","EAS",2166,197,0.974,0.989,True),("GTEx blood FCAR junction","EUR",3055,333,0.642,0.687,True)],
 columns=["trait","GWAS","variants_full","variants_truncated","PP.H4_full","PP.H4_truncated","same_top_variant"])
ws.cell(ws.max_row+2,1,'B. eQTL Catalogue targets').font = copy(hf); write_df(ws, B, row=ws.max_row+1)
rd = wb['README']; s = rd.max_row+1
for i,(a,b,c) in enumerate([
 ('S6 (updated)','Ancestry decomposition with bootstrap 95% CIs','Percentile intervals from 20 000 parametric bootstrap replicates; assumptions listed above the totals.'),
 ('S8 (updated)','Matched testing: exact versions','All genes; excluding immune-cell CFH gene tests; excluding the whole CFH window (22 windows).'),
 ('S9','Prior sensitivity (p12 from 1e-6 to 1e-4)','PP.H4 recomputed exactly; immune-cell window counts above p12 = 1e-5 are lower bounds.'),
 ('S10','Steiger directionality test','Variance explained by the instrument in the molecular trait vs IgAN.'),
 ('S12','Truncation benchmark','Complete ImmuNexUT and eQTL Catalogue data re-analyzed after keeping only variants with QTL P < 0.05.')]):
    rd.cell(s+i,1,a); rd.cell(s+i,2,b); rd.cell(s+i,3,c)
wb.save(sys.argv[3]); print(wb.sheetnames)
