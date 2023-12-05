import os
import gzip
import tarfile
import scipy.io
import numpy as np
import pandas as pd
import anndata as ad
import tacco as tc

tar_archive = snakemake.input['tar']
pelka_h5ad = snakemake.input['pelka']
output_h5ad = snakemake.output['h5ad']

with tarfile.open(tar_archive, mode="r") as tar:
    tarnames = tar.getnames()
    
    tarnames = pd.concat([pd.Series(tarnames).str.rsplit('_',n=1,expand=True).rename(columns={0:'samples',1:'extensions'}),pd.DataFrame({'tarnames':tarnames})],axis=1)
    
    sets = tarnames.pivot(index='samples', columns='extensions', values='tarnames')
    
    adatas = {}
    for sample in sets.index:
        with tar.extractfile(sets.loc[sample, 'barcodes.tsv.gz']) as f:
            with gzip.GzipFile(mode='r', fileobj=f) as gz:
                barcodes = pd.read_csv(gz, sep='\t', names=['bc'])
            barcodes['bc'] = barcodes['bc'].str.split('-1').str[0]
            barcodes = barcodes.set_index('bc')
            barcodes['sample'] = sample
        with tar.extractfile(sets.loc[sample, 'features.tsv.gz']) as f:
            with gzip.GzipFile(mode='r', fileobj=f) as gz:
                features = pd.read_csv(gz, sep='\t', names=['ENSG','GENE','feature_type'])
            features = features.set_index('ENSG')[['GENE']]
        with tar.extractfile(sets.loc[sample, 'matrix.mtx.gz']) as f:
            with gzip.GzipFile(mode='r', fileobj=f) as gz:
                matrix = scipy.io.mmread(gz).T.tocsr()
        adatas[sample] = ad.AnnData(matrix, obs=barcodes, var=features, dtype=np.float32)

adata = ad.concat(adatas.values(), index_unique='_', merge='same')
adata.var = adata.var.reset_index().set_index('GENE')
adata.var_names_make_unique()

for key,values in adata.obs['sample'].str.split('_',n=3,expand=True).rename(columns={0:'GSM',1:'Patient',2:'State',3:'ID'}).items():
    adata.obs[key] = values
adata.obs['ID'] = adata.obs['ID'].fillna('1')
adata = adata[adata.obs['State'] != 'blood'].copy()

adata.obs['sample'] = adata.obs['sample'].astype('category')

adata = tc.pp.filter(adata, min_genes_per_cell=200, )

# transfer epithelial compartment annotation from other dataset
human_pelka = ad.read(pelka_h5ad)
tc.tl.annotate(adata, human_pelka, 'cl295v11SubShort', result_key='pelka_cl295v11SubShort', verbose=0);
for derived_anno in ['compartment','clTopLevel', 'clMidwayPr', 'cl295v11SubFull',]:
    tc.utils.merge_annotation(adata,'pelka_cl295v11SubShort',mapping={k:v['cl295v11SubShort'].to_numpy() for k,v in human_pelka.obs[['cl295v11SubShort',derived_anno]].drop_duplicates().groupby(derived_anno)}, result_key=f'pelka_{derived_anno}');
adata.obs['pelka_epithelial'] = adata.obsm['pelka_compartment']['epithelial'] > 0.95

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')
