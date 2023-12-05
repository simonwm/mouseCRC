import os
import gzip
import tarfile
import scipy.io
import numpy as np
import pandas as pd
import anndata as ad
import tacco as tc

mtx_file = snakemake.input['matrix']
genes_file = snakemake.input['genes']
barcodes_file = snakemake.input['barcodes']
pelka_h5ad = snakemake.input['pelka']
output_h5ad = snakemake.output['h5ad']

X = scipy.io.mmread(mtx_file).T.tocsr()

features=pd.read_csv(genes_file, names=['ENSG','GENE'], sep='\t')
features = features.set_index('ENSG')[['GENE']]

barcodes=pd.read_csv(barcodes_file, names=['bc'], sep='\t')
barcodes = barcodes['bc'].str.split('_',expand=True).rename(columns={1:'patient',2:'State'})[['patient','State']].set_index(barcodes['bc'].to_numpy())
barcodes['sample'] = barcodes['patient'] + '_' + barcodes['State']
barcodes = barcodes[['sample','patient','State']]

adata = ad.AnnData(X, obs=barcodes, var=features, dtype=np.float32)

adata.var = adata.var.reset_index().set_index('GENE')
adata.var_names_make_unique()

adata = adata[adata.obs['State'].isin(['CRC', 'LM'])].copy()

adata.obs['sample'] = adata.obs['sample'].astype('category')

adata = tc.pp.filter(adata, min_genes_per_cell=500, )

# transfer epithelial compartment annotation from other datasets
human_pelka = ad.read(pelka_h5ad)
tc.tl.annotate(adata, human_pelka, 'cl295v11SubShort', result_key='pelka_cl295v11SubShort', verbose=0);
for derived_anno in ['compartment','clTopLevel', 'clMidwayPr', 'cl295v11SubFull',]:
    tc.utils.merge_annotation(adata,'pelka_cl295v11SubShort',mapping={k:v['cl295v11SubShort'].to_numpy() for k,v in human_pelka.obs[['cl295v11SubShort',derived_anno]].drop_duplicates().groupby(derived_anno)}, result_key=f'pelka_{derived_anno}');
adata.obs['pelka_epithelial'] = adata.obsm['pelka_compartment']['epithelial'] > 0.95

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')
