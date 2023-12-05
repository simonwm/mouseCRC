import os
import pandas as pd
import scanpy as sc
import tacco as tc

commot_sample_h5ads = snakemake.input

mouse_slideseq_pathway_h5ad = snakemake.output['mouse_slideseq_pathway_h5ad']

adatas = {}
samples = []
for commot_sample_h5ad in commot_sample_h5ads:
    sample = os.path.splitext(os.path.basename(commot_sample_h5ad))[0]
    adatas[sample] = sc.read(commot_sample_h5ad)
    samples.append(sample)

concated = sc.concat(adatas)

{sample:tc.testing.assert_frame_equal(adatas['AV10a_P'].uns['commot-cellchat-info']['df_ligrec'],adatas[sample].uns['commot-cellchat-info']['df_ligrec']) for sample in samples}
concated.uns['commot-cellchat-info'] = adatas['AV10a_P'].uns['commot-cellchat-info']

pathways = concated.uns['commot-cellchat-info']['df_ligrec']['pathway'].unique()

concated.obs['region'] = concated.obs['region'].astype('category')
concated.obs['State'] = concated.obs['State'].astype('category')
concated.obs['SampleID'] = concated.obs['SampleID'].astype('category')

rec_sums = concated.obsm['commot-cellchat-sum-receiver'][['r-'+p for p in pathways]]
sen_sums = concated.obsm['commot-cellchat-sum-sender'][['s-'+p for p in pathways]]
recsen_sums = pd.concat([rec_sums,sen_sums],axis=1)
adata_pathway = sc.AnnData(X=recsen_sums, obs=concated.obs)

os.makedirs(os.path.dirname(mouse_slideseq_pathway_h5ad), exist_ok=True)
adata_pathway.write(mouse_slideseq_pathway_h5ad, compression='gzip')
