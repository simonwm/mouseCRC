import os
import gzip
import tarfile
import scipy.io
import numpy as np
import pandas as pd
import scanpy as sc
import tacco as tc

input_h5 = snakemake.input['h5']
input_anno_cluster = snakemake.input['anno_cluster']
input_anno_meta = snakemake.input['anno_meta']
output_h5ad = snakemake.output['h5ad']
ccNMF_CellAct_dir = snakemake.input['ccNMF_CellAct']
ccNMF_GeneWeight_dir = snakemake.input['ccNMF_GeneWeight']

# load data
data = sc.read_10x_h5(input_h5)
cluster = pd.read_csv(input_anno_cluster).set_index('sampleID')
for c in cluster.columns:
    data.obs[c] = cluster[c]
metatables = pd.read_csv(input_anno_meta).set_index('cellID')
for c in metatables.columns:
    data.obs[c] = metatables[c]
# filter out non-genes
data = data[:,~(data.var.index.str.startswith('HASH-S') | data.var.index.str.startswith('ADT-'))].copy()

# create derived annotations
data.obs['TMMR'] = (data.obs['SPECIMEN_TYPE'].astype(str)+"-"+data.obs['MMRStatus'].astype(str)).map({'T-MMRp':'MMRp','T-MMRd':'MMRd','N-nan':'normal',}).astype('category')
data.obs['compartment'] = data.obs['clTopLevel'].map({'Epi':'epithelial','Strom':'stromal','TNKILC':'immune','Myeloid':'immune','B':'immune','Mast':'immune','Plasma':'immune'})
data.obs['PatientTypeID'] = data.obs['PatientTypeID'].astype('category')
data.obs['PID'] = data.obs['PID'].astype('category')
tc.utils.merge_annotation(data,'clTopLevel',{'epithelial':['Epi',],'stromal':['Strom',],'immune':['TNKILC','Myeloid','B','Mast','Plasma']},'compartment')

# load programs and weights
for prog in ['B','Epi','EpiTd','EpiTp','Mast','Myeloid','Plasma','Strom','T']:
    label = prog
    prog_name = f'ccNMF_{prog}'
    if prog == 'Epi':
        label = 'EpiTGlobalv5ForceK43'
    elif prog == 'EpiTd':
        label = 'EpiTMSIv4ForceK29'
    elif prog == 'EpiTp':
        label = 'EpiTMSSv4ForceK32'
    df = pd.read_csv(f'{ccNMF_CellAct_dir}/ccNMF_cell_{label}.csv.gz').set_index('Var1')
    data.obsm[prog_name] = df.reindex(index=data.obs.index)
    data.obsm[prog_name] /= data.obsm[prog_name].sum(axis=1).to_numpy()[:,None]
    df = pd.read_csv(f'{ccNMF_GeneWeight_dir}/ccNMF_RawWeights_{label}.csv.gz').set_index('Var1')
    data.varm[prog_name] = df.set_index(data.var.index)
    data.varm[prog_name] /= data.varm[prog_name].sum(axis=0).to_numpy()[None,:]

data.X = data.X.astype(np.float32)
data.var_names_make_unique()

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
data.write(output_h5ad, compression='gzip')
