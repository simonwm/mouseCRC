import os
import scanpy as sc
import commot as ct

mouse_slideseq_h5ad = snakemake.input['slideseq']

sample = snakemake.wildcards['sample']

commot_sample_h5ad = snakemake.output['commot_sample_h5ad']

tdata = sc.read(mouse_slideseq_h5ad)

def run_COMMOT(adata, x_key, y_key):
    adata = adata[adata.X.sum(axis=1)>=100].copy()
    
    adata.obsm['spatial'] = adata.obs[[x_key, y_key]].to_numpy()

    # Basic preprocessing as in https://commot.readthedocs.io/en/latest/notebooks/Basic_usage.html
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)

    # loading CellChat database as in https://commot.readthedocs.io/en/latest/notebooks/Basic_usage.html
    df_cellchat = ct.pp.ligand_receptor_database(database='CellChat', species='mouse')
    # filtering CellChat database as in slideseqv2-mouse-hippocampus/1-lr_signaling.ipynb from https://doi.org/10.5281/zenodo.7272562
    df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata, min_cell_pct=0.0)

    # reconstruct communication as in slideseqv2-mouse-hippocampus/1-lr_signaling.ipynb from https://doi.org/10.5281/zenodo.7272562
    ct.tl.spatial_communication(adata,
        database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=200, heteromeric=True, pathway_sum=True)
    
    return adata

cdata = run_COMMOT(tdata[tdata.obs['SampleID'].isin([sample])], x_key='x', y_key='y')

os.makedirs(os.path.dirname(commot_sample_h5ad), exist_ok=True)
cdata.write(commot_sample_h5ad, compression='gzip')
