import os
import gzip
import tarfile
import scipy.io
import numpy as np
import pandas as pd
import scanpy as sc
import tacco as tc

gex_file = snakemake.input['gex']
anno_file = snakemake.input['anno']
coadread_anno_file = snakemake.input['coadread_anno']
output_h5ad = snakemake.output['h5ad']

TCGA_X = pd.read_csv(gex_file, sep='\t', index_col=0).T
TCGA_obs = pd.read_excel(anno_file, sheet_name=0, index_col=0)
TCGA_coadread = pd.read_csv(coadread_anno_file, sep='\t', index_col=0)

adata = sc.AnnData(TCGA_X, obs=TCGA_obs.set_index(TCGA_obs['bcr_patient_barcode']).reindex(TCGA_X.index.str.rsplit('-',n=4).str[0]).set_index(TCGA_X.index), dtype=np.float32)

TCGA_coadread['MSI'] = TCGA_coadread['CDE_ID_3226963']
still_na = TCGA_coadread['MSI'].isna()
TCGA_coadread.loc[still_na,'MSI'] = TCGA_coadread.loc[still_na,'MSI_updated_Oct62011']
changing = (~TCGA_coadread['CDE_ID_3226963'].isna()) & (~TCGA_coadread['MSI_updated_Oct62011'].isna()) & (TCGA_coadread['CDE_ID_3226963'] != TCGA_coadread['MSI_updated_Oct62011'])
TCGA_coadread.loc[changing,'MSI'] = 'Indeterminate'
patient_MSIstatus_map = TCGA_coadread[['MSI','_PATIENT']].drop_duplicates().set_index('_PATIENT')['MSI']
adata.obs['MSIstatus'] = adata.obs['bcr_patient_barcode'].map(patient_MSIstatus_map)

TCGA_genes = adata.var_names.str.split('|').str[0]
multiplicity = TCGA_genes.value_counts()
duplicates = multiplicity[multiplicity>1].index
adata = adata[:,~TCGA_genes.isin(duplicates)].copy()
adata.var.index = adata.var_names.str.split('|').str[0]

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')
