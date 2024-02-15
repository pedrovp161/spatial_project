import scanpy as sc
import pandas as pd
from anndata import AnnData
from PIL import Image
import numpy as np
import json
import os

def read_free_h5ad(DIR):
    # File paths
    pos_path = os.path.join(DIR, "spatial", "tissue_positions_list.csv")
    bar_path = os.path.join(DIR, "raw_feature_bc_matrix", "barcodes.tsv.gz")
    feat_path = os.path.join(DIR, "raw_feature_bc_matrix", "features.tsv.gz")
    matx_path = os.path.join(DIR, "raw_feature_bc_matrix", "matrix.mtx.gz")
    json_path = os.path.join(DIR, "spatial", "scalefactors_json.json")

    # Finding image files
    spatial_path = os.path.join(DIR, "spatial")
    hier_files = [file for file in os.listdir(spatial_path) if file.endswith("hires_image.png")]
    lower_files = [file for file in os.listdir(spatial_path) if file.endswith("lowres_image.png")]
    hier_path = os.path.join(spatial_path, hier_files[0])  # Assuming there's only one hires image
    lower_path = os.path.join(spatial_path, lower_files[0])  # Assuming there's only one lowres image

    # transforming image into numpy array
    imh = Image.open(hier_path)
    iml = Image.open(lower_path)
    image_hirer = np.array(imh)
    image_lower = np.array(iml)

    # Reading the matrix file and transposing it
    adata = sc.read(matx_path, delimiter='\t')
    adata = adata.T

    # Setting up gene information in adata.var
    genes = pd.read_csv(feat_path, delimiter="\t", header=None, index_col=0)
    genes.index.name = None
    genes = genes.rename(columns={1:'gene_ids', 2:'feature_types'})
    adata.var = genes

    # Reading the barcodes
    barcodes = pd.read_csv(bar_path, header=None, delimiter="\t")
    barcodes.index.name = None

    # Reading positional information
    pos = pd.read_csv(pos_path, header=None)
    
    # Merging positional information with barcodes
    pos = pd.merge(pos, barcodes, how='inner', left_on=0, right_on=0)
    
    # Selecting and renaming columns for positional information
    pos = pos[[0, 4, 5]]
    pos.index = pos[0]
    del pos[0]
    pos = pos.rename(columns={4:"array_row", 5:"array_col"})
    pos.index.name = None

    # Setting up the adata.obs
    adata.obs = pos

    # Setting up the adata.obsm["spatial"]
    coordinates = pos.values

    # Making the uns model
    with open(json_path, "r") as arquivo_json:
        scale_info = json.load(arquivo_json)
    modelo_uns = {
        'spatial': {
            'library_id': {
                'images': {
                    'hires': image_hirer,
                    'lowres': image_lower
                },
                'scalefactors': scale_info,
                'metadata': {
                    'chemistry_description': "Spatial 3' v1",
                    'software_version': 'spaceranger-1.2.0'
                }
            }
        }
    }

    adata.uns = modelo_uns  # Assigning the uns model to adata.uns

    # spatial coordinates
    adata.obsm["spatial"] = pos.values

    return adata



