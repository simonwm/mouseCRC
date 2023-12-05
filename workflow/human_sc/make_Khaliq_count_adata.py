import os
import numpy as np
import scipy.sparse
import pandas as pd
import anndata as ad
import tacco as tc

counts_file = snakemake.input['counts']
anno_file = snakemake.input['anno']
pelka_h5ad = snakemake.input['pelka']
output_h5ad = snakemake.output['h5ad']

df = pd.read_csv(counts_file, index_col=0)
df_anno = pd.read_csv(anno_file, index_col=0)
adata = ad.AnnData(df.T, obs=df_anno, dtype=np.float32)
adata.X = scipy.sparse.csr_matrix(adata.X)

# transfer epithelial compartment annotation from other dataset
human_pelka = ad.read(pelka_h5ad)
tc.tl.annotate(adata, human_pelka, 'cl295v11SubShort', result_key='pelka_cl295v11SubShort', verbose=0);
for derived_anno in ['compartment','clTopLevel', 'clMidwayPr', 'cl295v11SubFull',]:
    tc.utils.merge_annotation(adata,'pelka_cl295v11SubShort',mapping={k:v['cl295v11SubShort'].to_numpy() for k,v in human_pelka.obs[['cl295v11SubShort',derived_anno]].drop_duplicates().groupby(derived_anno)}, result_key=f'pelka_{derived_anno}');
adata.obs['pelka_epithelial'] = adata.obsm['pelka_compartment']['epithelial'] > 0.95

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')
