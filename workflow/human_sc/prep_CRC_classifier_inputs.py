# Prepare normalized pseudo-bulk data for human and mouse datasets for use in CRC classification

import os
import warnings
warnings.filterwarnings('ignore','invalid value encountered in true_divide')

import pandas as pd
import numpy as np
import anndata as ad
import scipy

import sys
# Make helper functions available: The notebook expects to be executed either in the sub-workflow directory or in the notebooks directory
#sys.path.insert(1, '../../'), sys.path.insert(1, '../'); # prefer to look just one directory up
sys.path.insert(1, 'workflow/')
import helper
sys.path.pop(1), sys.path.pop(1);

get_path = helper.get_paths('human_sc')

import tacco as tc

CMS_classifier_input = snakemake.output['CMS_classifier_input']

## Load data and convert the genes to the common MGI homology classes

tc.tl.setup_orthology_converter(f'{get_path("resources")}/MGI/HOM_AllOrganism.rpt');

mouse_scrna = ad.read(f'{get_path("resources","mouse_sc")}/scRNAseq.h5ad')
mouse_slideseq = ad.read(f'{get_path("resources","mouse_slideseq")}/slideseq.h5ad')
mouse_slideseq_by_compartment = ad.read(f'{get_path("resources","mouse_slideseq")}/slideseq_by_compartment.h5ad')

# use only data from sufficiently covered beads
mouse_slideseq = mouse_slideseq[mouse_slideseq.X.sum(axis=1)>=100].copy()
mouse_slideseq_by_compartment = mouse_slideseq_by_compartment[mouse_slideseq_by_compartment.obs['index'].isin(mouse_slideseq.obs.index)].copy()

mouse_scrna = tc.tl.run_orthology_converter(mouse_scrna, 'mouse', target_tax_id='human', use_synonyms=True)
mouse_slideseq_by_compartment = tc.tl.run_orthology_converter(mouse_slideseq_by_compartment, 'mouse', target_tax_id='human', use_synonyms=True)
# remove varm as it breaks anndata if subsetted...
for k in [k for k in mouse_scrna.varm]:
    del mouse_scrna.varm[k]

mouse_adatas = { 'scRNA': mouse_scrna, 'SlideSeq': mouse_slideseq_by_compartment }

human_data_sources = [ 'Pelka', 'Becker', 'Zheng', 'Khaliq', 'Che', 'Chen', 'Joanito_3p', 'Joanito_5p' ]

human_adatas = { source: ad.read(f'{get_path("data")}/{source}.h5ad') for source in human_data_sources }

for source, human_adata in human_adatas.items():
    human_adata.obs['species'] = 'human'
    human_adata.obs['source'] = source

for source, mouse_adata in mouse_adatas.items():
    mouse_adata.obs['species'] = 'mouse'
    mouse_adata.obs['source'] = source

all_adatas = {**mouse_adatas, **human_adatas}

# generate a consistent layer of annotation across datasets
for source, adata in all_adatas.items():
    if 'SampleID' not in adata.obs:
        if 'samples' in adata.obs:
            adata.obs['SampleID'] = adata.obs['samples']
        elif 'sample' in adata.obs:
            adata.obs['SampleID'] = adata.obs['sample']
        elif 'PID' in adata.obs:
            adata.obs['SampleID'] = adata.obs['PID']
        else:
            raise ValueError(f'Unknown sample spec in source {source}!')
    if 'State' not in adata.obs:
        if 'GrossPathology' in adata.obs:
            adata.obs['State'] = adata.obs['GrossPathology']
        elif 'TMMR' in adata.obs:
            adata.obs['State'] = adata.obs['TMMR']
        elif 'Condition' in adata.obs:
            adata.obs['State'] = adata.obs['Condition']
        else:
            raise ValueError(f'Unknown state spec in source {source}!')
    if 'Epithelial' not in adata.obs:
        if 'epithelial' in adata.obs:
            adata.obs['Epithelial'] = adata.obs['epithelial']
        elif 'compartment' in adata.obs:
            adata.obs['Epithelial'] = adata.obs['compartment'].isin(['Epi','Epithelial','epi','epithelial'])
        elif 'pelka_epithelial' in adata.obs:
            adata.obs['Epithelial'] = adata.obs['pelka_epithelial']
        else:
            raise ValueError(f'Unknown sample spec in source {source}!')
    adata.obs['source-State'] = adata.obs['source'].astype(str) + ': ' + adata.obs['State'].astype(str)

# filter all datasets to a common gene set
filtered_adatas = tc.pp.filter([*human_adatas.values(),*mouse_adatas.values()], return_view=False, remove_constant_genes=True)

def prep_pseudobulk(adata):
    pseudobulk = tc.tl.merge_observations(adata, 'SampleID')
    pseudobulk.obs['SampleID'] = pseudobulk.obs['SampleID'].astype(str)
    pseudobulk.obs = pseudobulk.obs.set_index('SampleID')
    return pseudobulk

# generate full pseudobulk profiles for the CMS classifier

all_pseudobulk = [ prep_pseudobulk(adata) for adata in filtered_adatas ]

adata_pseudobulk = ad.concat(all_pseudobulk)

# normalize data for CMS classifier

pseudobulk = adata_pseudobulk.to_df().T # transpose for R style

pseudobulk = 1e6 * pseudobulk / pseudobulk.sum(axis=1).to_numpy()[:,None]

pseudobulk = np.log2(pseudobulk+1)

# translate and export data for CMS classifier

HOM = pd.read_csv(f'{get_path("resources")}/MGI/HOM_AllOrganism.rpt',sep='\t')

gene_name2EntrezID = HOM[HOM['Common Organism Name']=='human'].set_index('Symbol')['EntrezGene ID']

mapping_multiplicities = gene_name2EntrezID.index.value_counts()
gene_name2EntrezID = gene_name2EntrezID[mapping_multiplicities[mapping_multiplicities==1].index]

pseudobulk = pseudobulk.loc[pseudobulk.index.isin(gene_name2EntrezID.index),:]

pseudobulk.index = pseudobulk.index.map(gene_name2EntrezID)

pseudobulk.to_csv(CMS_classifier_input)

