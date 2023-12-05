import os
import gzip
import tarfile
import scipy.io
import numpy as np
import pandas as pd
import anndata as ad
import tacco as tc
import tables

epi_count_h5 = snakemake.input['epi_count']
epi_meta_csv = snakemake.input['epi_meta']
nonepi_count_h5 = snakemake.input['nonepi_count']
nonepi_meta_csv = snakemake.input['nonepi_meta']
clinical_meta_csv = snakemake.input['clinical_meta']

output_h5ad_3p = snakemake.output['h5ad_3p']
output_h5ad_5p = snakemake.output['h5ad_5p']

def read_h5(filename):
    with tables.open_file(filename, "r") as h5file:
        
        barcodes = [s.decode('ascii') for s in h5file.root.matrix.barcodes.read()]
        data = h5file.root.matrix.data.read()
        indices = h5file.root.matrix.indices.read()
        indptr = h5file.root.matrix.indptr.read()
        shape = h5file.root.matrix.shape.read()
        name = [s.decode('ascii') for s in h5file.root.matrix.features.name.read()]
        
    X = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
    
    adata = ad.AnnData(X.T, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=name), dtype=np.float32)
    
    return adata

epi_count = read_h5(epi_count_h5)
nonepi_count = read_h5(nonepi_count_h5)
epi_meta = pd.read_csv(epi_meta_csv).set_index('cell.ID')
nonepi_meta = pd.read_csv(nonepi_meta_csv).set_index('cell.ID')

nonepi_meta['iCMS'] = np.nan
nonepi_meta['msi'] = nonepi_meta['sample.ID'].map(epi_meta[['sample.ID','msi']].drop_duplicates().set_index('sample.ID')['msi'])

for c in epi_meta.columns:
    epi_count.obs[c] = epi_meta[c]

for c in nonepi_meta.columns:
    nonepi_count.obs[c] = nonepi_meta[c]

adata = ad.concat({True:epi_count,False:nonepi_count},label='Epithelial')
adata.obs['Epithelial'] = adata.obs['Epithelial'].astype(bool)

adata = adata[~adata.obs['sample.origin'].isin(['LymphNode'])].copy()

adata.obs['State'] = adata.obs['sample.origin']
adata.obs.loc[(adata.obs['sample.origin'].str.startswith('Tumor')),'State'] = adata.obs['msi']
adata.obs['SampleID'] = adata.obs['sample.ID']
adata.obs['SampleID'] = adata.obs['SampleID'].astype('category')

adata_3p = adata[adata.obs['dataset'].isin(['CRC-SG2', 'SMC', 'KUL3',])].copy()
adata_5p = adata[adata.obs['dataset'].isin(['CRC-SG1', 'KUL5',])].copy()

os.makedirs(os.path.dirname(output_h5ad_3p), exist_ok=True)
adata_3p.write(output_h5ad_3p, compression='gzip')

os.makedirs(os.path.dirname(output_h5ad_5p), exist_ok=True)
adata_5p.write(output_h5ad_5p, compression='gzip')