"""
    Utility functions for use in notebooks related to prepping the modality quantifications
"""

# imports
from anndata import AnnData
from pandas import DataFrame
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import json
import numpy as np

# constants


# data quick peek functions
def peek_anndata(adata: AnnData, message: str=None, verbose: bool=False):
    if not message is None and len(message) > 0:
        print(message)
    print(adata)
    if verbose:
        display(adata.obs.head())
        display(adata.var.head())

def peek_dataframe(df: DataFrame, verbose: bool=False):
    print(f'{df.shape=}')
    if verbose:
        display(df.head())

def scanpy_umap(this_adata: AnnData, this_feature: str, proj_name: str, out_dir: str=None, 
                legend_on: bool=True, this_dpi: int=100, width: int=9, height: int=9):
    # if out_dir specified then format figure file name
    if not out_dir is None and len(out_dir) > 0:
        sc.settings.figdir = f'{out_dir}/'
        figure_file = f'_{proj_name}.umap.{this_feature}.png'
    else:
        figure_file = None
    # setup if legend is on or off the figure
    if legend_on:
        legend_param = 'on data'
    else:
        legend_param = 'right margin'
    # generate the figures
    with rc_context({'figure.figsize': (width, height), 'figure.dpi': this_dpi}):
        plt.style.use('seaborn-v0_8-talk')
        sc.pl.umap(this_adata, color=[this_feature], 
                   frameon=False, legend_loc=legend_param, save=figure_file)

def scanpy_dotplot(this_adata: AnnData, this_feature: str, proj_name: str, 
                   this_dict: dict, out_dir: str=None, layer_name: str=None,
                   this_dpi: int=100, width: int=9, height: int=9):
    # if out_dir specified then format figure file name
    if not out_dir is None and len(out_dir) > 0:
        sc.settings.figdir = f'{out_dir}/'
        figure_file = f'_{proj_name}.markers_dotplot.{this_feature}.png'
    else:
        figure_file = None
    # generate the figures
    with rc_context({'figure.figsize': (width, height), 'figure.dpi': this_dpi}):
        plt.style.use('seaborn-v0_8-talk')
        sc.pl.dotplot(this_adata, this_dict, groupby=this_feature, cmap='Purples', 
                      layer=layer_name, save=figure_file)        

def load_cell_type_markers(marker_genes_json: str, adata: AnnData, 
                          filter_hv: bool=True, verbose: bool=False) -> (set, dict):
    markers_dict = None
    markers = None
    if filter_hv:
        # get the high variance feature set
        possible_features = set(adata.var.loc[adata.var.highly_variable].index.values)
        print(f'{len(possible_features)} features are labeled high variance')
    else:
        possible_features = set(adata.var.index.values)
        print(f'{len(possible_features)} features are present')        
    
    with open(marker_genes_json, 'r') as json_file:
        markers_dict = json.load(json_file)
    # get the set of all markers across the the cell-types
    markers = {item for sublist in markers_dict.values() for item in sublist}
    print(f'{len(markers)} marker features loaded')
    # find the marker genes that are present in the current HV features
    missing_markers = markers - possible_features
    print(f'missing {len(missing_markers)} markers: {missing_markers}')
    # drop the markers missing for the current HV features
    markers = markers & possible_features
    print(f'{len(markers)} marker features found')
    if verbose:
        print(f'markers found: {markers}')
    
    # update cell-type markers dict to drop any of the missing markers
    list_keys_to_delete = []
    for cell_type, marker_list in markers_dict.items():
        new_list = list(set(marker_list) & markers)
        if len(new_list) > 0:
            markers_dict[cell_type] = new_list
        else:
            list_keys_to_delete.append(cell_type)
    for cell_type in list_keys_to_delete:
        markers_dict.pop(cell_type)
    
    return markers, markers_dict

def heatmap_compare(adata: AnnData, set1: str, set2: str):
    this_df = (
        adata.obs.groupby([set1, set2])
        .size()
        .unstack(fill_value=0)
    )
    norm_df = this_df/this_df.sum(axis=0)

    with rc_context({'figure.figsize': (12, 12), 'figure.dpi': 100}):
        plt.style.use('seaborn-v0_8-bright')
        _ = plt.pcolor(norm_df, edgecolor='black')
        _ = plt.xticks(np.arange(0.5, len(this_df.columns), 1), this_df.columns, rotation=90)
        _ = plt.yticks(np.arange(0.5, len(this_df.index), 1), this_df.index)
        plt.xlabel(set2)
        plt.ylabel(set1)
        plt.show()