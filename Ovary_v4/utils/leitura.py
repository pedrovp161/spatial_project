import scanpy as sc
import numpy as np
import pandas as pd
import os
from PIL import Image
import json


class Leitura():
    def __init__(self, dir) -> None:
        # definindo diretorio
        self.DIR = dir
        self.pasta_raw_GEO = self.DIR + "/data/GEO_total"
        self.pasta_filtrado_GEO = self.DIR + "/data/Resultados/GEO_files"

    def listar_raw_GEO_dic_com_dados(self):# Interface para dados raw
        return self.listar_dic_com_dados(self.pasta_raw_GEO)

    def listar_filtrado_GEO_dic_com_dados(self):# Interface para criar dic com pasta de dados filtrados
        return self.listar_dic_com_dados(self.pasta_filtrado_GEO)

    def listar_caminhos_arquivos(self, pasta):

        caminhos_arquivos = []

        # Lista os nomes de arquivos na pasta
        nomes_arquivos = os.listdir(pasta)

        # Cria os caminhos completos para cada arquivo
        for nome_arquivo in nomes_arquivos:
            caminho_completo = os.path.join(pasta, nome_arquivo)
            caminhos_arquivos.append(caminho_completo)

        return caminhos_arquivos
    
    def listar_dic_com_dados(self, pasta):
        # Cria dicionários para armazenar os resultados
        paths = {}
        
        caminhos_subpastas = [os.path.join(pasta, nome) for nome in os.listdir(pasta) if os.path.isdir(os.path.join(pasta, nome))]

        # Itera sobre as subpastas e aplica a função read_free_h5ad ou sc.read_h5ad
        for caminho_subpasta in caminhos_subpastas:
            nome_subpasta = os.path.basename(caminho_subpasta)

            print("Processando subpasta:", nome_subpasta)  # Adicionando verificação de diretório

            if pasta == self.pasta_raw_GEO:
                paths[nome_subpasta] = self.read_free_h5ad(caminho_subpasta)
            elif pasta == self.pasta_filtrado_GEO:
                paths[nome_subpasta] = sc.read_h5ad(caminho_subpasta)
            else:
                print("ERROR")
        
        return paths
    
    def read_free_h5ad(self, DIR):
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
    



if __name__ == "__main__":# teste da classe
    import pathlib

    DIR = pathlib.Path(os.path.dirname(__file__)).parent #parent utilizado aqui para funcionar como teste da classe

    ler = Leitura(str(DIR))
