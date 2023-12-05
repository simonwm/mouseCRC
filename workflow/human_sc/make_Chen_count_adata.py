import os
import gzip
import tarfile
import scipy.io
import numpy as np
import pandas as pd
import anndata as ad
import tacco as tc

abnormal_epi_h5ad = snakemake.input['abnormal_epi_h5ad']
dis_epi_h5ad = snakemake.input['dis_epi_h5ad']
valdis_nonepi_h5ad = snakemake.input['valdis_nonepi_h5ad']
val_epi_h5ad = snakemake.input['val_epi_h5ad']
output_h5ad = snakemake.output['h5ad']

abnormals_epi = ad.read(abnormal_epi_h5ad)
dis_epi = ad.read(dis_epi_h5ad)
val_dis_nonepi = ad.read(valdis_nonepi_h5ad)
val_epi = ad.read(val_epi_h5ad)

abnormals_epi.X = scipy.sparse.csr_matrix(abnormals_epi.X)
dis_epi.X = scipy.sparse.csr_matrix(dis_epi.X)
val_dis_nonepi.X = scipy.sparse.csr_matrix(val_dis_nonepi.X)
val_epi.X = scipy.sparse.csr_matrix(val_epi.X)

abnormals_epi.obs['Cell_Type'] = 'abnormal'
abnormals_epi.obs['Epithelial'] = True
dis_epi.obs['Epithelial'] = True
val_epi.obs['Epithelial'] = True
val_dis_nonepi.obs['Epithelial'] = False

adata = ad.concat([abnormals_epi,dis_epi,val_dis_nonepi,val_epi])

adata.obs['State'] = adata.obs['Sample_Classification']
adata.obs['sample'] = adata.obs['HTAN Specimen ID']

adata.obs['sample'] = adata.obs['sample'].astype('category')

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')
