# helper functions for all sub-workflows

import os
import pandas as pd
import tacco as tc
import scipy
import numpy as np
import anndata as ad
import matplotlib

# general visualization settings
def apply_visualization_settings():
    highres = True
    default_dpi = 72.0 # matplotlib.rcParams['figure.dpi']
    if highres:
        matplotlib.rcParams['figure.dpi'] = 8 * default_dpi
    else:
        matplotlib.rcParams['figure.dpi'] = 2 * default_dpi
apply_visualization_settings()

def get_paths(subworkflow):
    # The code expects to be executed in the sub-workflow directory, the notebooks directory, or the base directory
    results_path = '../results'
    if not os.path.exists(results_path):
        results_path = f'../{results_path}'
    if not os.path.exists(results_path):
        results_path = f'results'
    
    if not os.path.exists(results_path):
        raise ValueError(f'The results path cannot be found!')
    
    def gen_paths(subworkflow):
        paths = {
            'results': results_path,
            'data': f'{results_path}/{subworkflow}/data',
            'plots': f'{results_path}/{subworkflow}/plots',
            'resources': f'{results_path}/../resources/{subworkflow}',
        }
        return paths
    default_paths = gen_paths(subworkflow)
    
    os.makedirs(default_paths['plots'], exist_ok=True)
    
    def get_path(which=None, subworkflow=None):
        if which == None:
            which = 'results'
        
        if subworkflow == None:
            paths = default_paths
        else:
            paths = gen_paths(subworkflow)
            
        if which in paths:
            return paths[which]
        else:
            raise ValueError(f'The path {which!r} is unkown!')
    
    return get_path

def get_colors(which):
    if which == 'mostly_gray_region':
        def make_mostly_gray_region_colors(region_colors):
            not_malignant = [t for t,c in region_colors.items() if 'Malignant-like' not in t]
            vals = np.linspace(0.3,0.9,len(not_malignant))[::-1]
            grays = vals.copy()
            #grays[0::2] = vals[len(vals)//2:]
            grays[0::2] = [vals[len(vals)//2],vals[len(vals)//2:][-1],*vals[len(vals)//2:][1:-1]]
            grays[1::2] = vals[:len(vals)//2]
            not_malignant_colors = {t:c for t,c in zip(not_malignant,[[v,v,v] for v in grays])}
            malignant_colors = {t:c for t,c in region_colors.items() if 'Malignant-like' in t }
            mostly_gray_region_colors = {**not_malignant_colors,**malignant_colors}
            return mostly_gray_region_colors
        return make_mostly_gray_region_colors(get_colors('region'))
    return {
        'compartment': tc.pl.get_default_colors(['epithelial','stromal','immune',]),
        'labels': {'Epi': (0.00784313725490196, 0.24313725490196078, 1.0), 'B': (0.10196078431372549, 0.788235294117647, 0.2196078431372549), 'TNK': (1.0, 0.48627450980392156, 0.0), 'Mono': (0.5490196078431373, 0.03137254901960784, 0.0), 'Mac': (0.9098039215686274, 0.0, 0.043137254901960784), 'Gran': (0.34901960784313724, 0.11764705882352941, 0.44313725490196076), 'Mast': (0.23529411764705882, 0.23529411764705882, 0.23529411764705882), 'Endo': (0.8549019607843137, 0.5450980392156862, 0.7647058823529411), 'Fibro': (0.6235294117647059, 0.2823529411764706, 0.0)},
        'cluster': {'Epi01 (Dysplastic Stem Like)': '#3c3c3c', 'Epi02 (Enterocytes)': '#592f0d', 'Epi03 (Stem/Progenitors)': '#b8850a', 'Epi04 (Secretory)': '#ff9f9b', 'Epi05 (Dysplastic Secretory Like)': '#d0bbff', 'Epi06': '#006374', 'Epi07 (Enteroendocrine)': '#a23582', 'Epi08': '#a1c9f4', 'Epi09': '#ffb482', 'Epi10': '#8de5a1', 'B01': '#023eff', 'B02': '#ff7c00', 'TNK01 (GdT/Cd8)': '#8c8c8c', 'TNK02 (Th1/Th17)': '#ccb974', 'TNK03 (NaiveTC)': '#64b5cd', 'TNK04 (GdT/Cd8)': '#001c7f', 'TNK05 (GdT_Il17+)': '#b1400d', 'TNK06 (Tregs)': '#12711c', 'TNK07 (NK)': '#da8bc3', 'TNK08 (ProliferatingTC)': '#8c0800', 'TNK09 (CD4_Tgfb1+)': '#591e71', 'Mono01': '#55a868', 'Mono02 (Dysplasia-Associated)': '#c44e52', 'Mono03 (Dysplasia-Associated, IFN)': '#8172b3', 'Mono04 (DCs)': '#937860', 'Mac01': '#00d7ff', 'Mac02 (Lyve1+)': '#4c72b0', 'Gran01': '#a3a3a3', 'Gran02': '#ffc400', 'Mast01': '#dd8452', 'Endo01 (Vascular)': '#1ac938', 'Endo02 (Lymphatic)': '#e8000b', 'Fibro01': '#8b2be2', 'Fibro02 (Myofibroblasts)': '#9f4800', 'Fibro03': '#f14cc1'},
        'program': {'Program 0': (0.00784313725490196, 0.24313725490196078, 1.0), 'Program 1 (Innate immune response)': (1.0, 0.48627450980392156, 0.0), 'Program 2 (Enteroendocrine)': (0.10196078431372549,  0.788235294117647,  0.2196078431372549), 'Program 3 (Proliferation (G2/M))': (0.9098039215686274,  0.0,  0.043137254901960784), 'Program 4 (Wnt signaling)': (0.5450980392156862,  0.16862745098039217,  0.8862745098039215), 'Program 5 (Transmembrane transport/Basolateral plasma membrane)': (0.6235294117647059,  0.2823529411764706,  0.0), 'Program 6 (Inflammatory response)': (0.9450980392156862,  0.2980392156862745,  0.7568627450980392), 'Program 7 (Innate immune response)': (0.6392156862745098,  0.6392156862745098,  0.6392156862745098), 'Program 8 (Apical plasma membrane)': (1.0, 0.7686274509803922, 0.0), 'Program 9': (0.0, 0.8431372549019608, 1.0), 'Program 10 (Oxidation-reduction process)': (0.0,  0.10980392156862745,  0.4980392156862745), 'Program 11 (Proliferation (G1/S))': (0.6941176470588235,  0.25098039215686274,  0.050980392156862744), 'Program 12': (0.07058823529411765, 0.44313725490196076, 0.10980392156862745), 'Program 13 (Deep crypt)': (0.5490196078431373, 0.03137254901960784, 0.0), 'Program 14 (Angiogenesis)': (0.34901960784313724,  0.11764705882352941,  0.44313725490196076), 'Program 15 (Goblet)': (0.34901960784313724,  0.1843137254901961,  0.050980392156862744), 'Program 16 (Stem cells)': (0.6352941176470588,  0.20784313725490197,  0.5098039215686274), 'Program 17': (0.23529411764705882, 0.23529411764705882, 0.23529411764705882), 'Program 18 (MHC II)': (0.7215686274509804,  0.5215686274509804,  0.0392156862745098), 'Program 19': (0.0, 0.38823529411764707, 0.4549019607843137)},
        'region': {'Region 0': '#8172b3', 'Region 1 (Stem cell niche - dysplastic)': '#8b2be2', 'Region 2 (Muscularis)': '#9f4800', 'Region 3 (Normal)': '#1ac938', 'Region 4': '#ffc400', 'Region 5 (Stem cell niche - normal)': '#f14cc1', 'Region 6 (Malignant-like, inflammation)': '#ff7c00', 'Region 7 ': '#00d7ff', 'Region 8 (Malignant-like, Deep crypt-like)': '#e8000b', 'Region 9 ': '#4c72b0', 'Region 10 (Normal)': '#dd8452', 'Region 11 (Malignant-like, EMT)': '#023eff', 'Region 12 (Normal)': '#55a868', 'Region 13': '#c44e52'},
    }[which]

def marker_genes(adata, group_key, rungo=True, reference_group=None, restrict_groups=None, goa_working_directory='.'):
    p_corr='fdr_bh'
    method='fisher'
    enrichments = tc.tl.enrichments(
        adata=adata,
        value_key=None,
        group_key=group_key,
        reference_group=reference_group,
        restrict_groups=restrict_groups,
        value_location='X',
        reduction=None,
        normalization=None,
        fillna=None,
        reads=False,
        min_obs= 1,
        p_corr=p_corr,
        method=method,
        direction='enrichment',
    )
    p_key = f'p_{method}_{p_corr}'
    enrichments = enrichments[(enrichments[p_key]<0.05)&(enrichments['enrichment']=='enriched')]
    
    if not rungo:
        return enrichments

    gene_label = 'value' if adata.var.index.name is None else adata.var.index.name
    if len(enrichments) == 0:
        print(f'no significant gene enrichments')
        return enrichments, None
    else:
        gene_lists = {region:list(df.sort_values(p_key)[gene_label]) for region, df in enrichments.groupby(group_key) if len(df)>0}
        gene_lists = { region: genes for region,genes in gene_lists.items() if len(genes) > 0}
        if len(gene_lists) > 0:
            tc.tl.setup_goa_analysis(adata.var.index, working_directory=goa_working_directory);
            results = tc.tl.run_goa_analysis(gene_lists)
            return enrichments, results
        else:
            print(f'no significant GO term enrichments')
            return enrichments, None

def map_short(x,i=0, split=' ', join=' '):
    if isinstance(x,str):
        return join.join(x.split(split)[:(i+1)])
    else:
        return {map_short(k,i,split=split,join=join):v for k,v in x.items()}

def annotate_epithelial_domain(adata):
    tc.tl.annotation_coordinate(adata, 'compartment', 'sample', result_key='comp_dist',max_distance=100,delta_distance=10, critical_neighbourhood_size=4.0,sparse=True, verbose=0);
    adata.obs['epi_domain'] = adata.obsm['comp_dist']['stromal'] > 75

def map_genes_for_cartana(adata):
    Ly6c1 = adata.var.index=='Ly6c1'
    Ly6c2 = adata.var.index=='Ly6c2'
    Trgc2 = adata.var.index=='Tcrg-C2'
    Cd244 = adata.var.index=='Cd244'
    renamed = adata.var.index.to_numpy()
    renamed[Ly6c1] = 'Ly6c1/2'
    renamed[Trgc2] = 'Trgc2'
    renamed[Cd244] = 'Cd244a'
    for k in adata.varm:
        adata.varm[k][Ly6c1] += adata.varm[k][Ly6c2]
        adata.varm[k].set_index(renamed,inplace=True)
    adata.X[Ly6c1] += adata.X[Ly6c2]
    adata.var.set_index(renamed,inplace=True)
    adata = adata[:,~Ly6c2].copy()
    return adata
