import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from PIL import Image

class Analysis:
    def __init__(self):
        self.data_dir = "C:\\Users\\pedro\\OneDrive\\Área de Trabalho\\projeto_INCA\\spatial_ovary\\Ovary_v4\\data"
        self.tenx_raw_dir = self.data_dir+"\\10x2"
        self.tenx_filtered_dir = self.data_dir+"\\10x2_filtered.h5ad"
        self.result_dir = self.data_dir+"\\Resultados"
        self.GEO_files = self.result_dir+"\\GEO_files"
        self.adatas = None
        self.adata_concat_GEO = None
        self.tenx_adata = None
        self.filenames = None
        
    def Filenames(self):

        self.filenames = os.listdir(self.GEO_files)
        order = [0,4,5,6,7,8,9,10,11,1,2,3]
        self.filenames = [self.filenames[i] for i in order]
        return self.filenames
    
    def load_all(self):
        self.adatas = [sc.read(str(self.GEO_files)+ "/" + i) for i in self.filenames]
        self.adata_concat_GEO = self.adatas[0].concatenate(self.adatas[1:], join = 'outer', fill_value = 0, uns_merge="first")
        self.tenx_adata = sc.read_h5ad(self.tenx_filtered_dir)
        return self.tenx_adata, self.adata_concat_GEO
    
    def plot_tx(self):
        self.Filenames()
        self.load_all()
        print("escreva qual dado vc deseja acessar ou help para instruções")
        pedido = input("imagem desejada: ")
        match pedido.lower:
            case "mt":
                sc.pl.spatial(self.tenx_adata, color=["pct_counts_mt"])
            case "graf":
                fig, axs = plt.subplots(1, 4, figsize=(15, 4))
                sns.distplot(self.tenx_adata.obs["total_counts"], kde=False, ax=axs[0])
                sns.distplot(self.tenx_adata.obs["total_counts"][self.tenx_adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
                sns.distplot(self.tenx_adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
                sns.distplot(self.tenx_adata.obs["n_genes_by_counts"][self.tenx_adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
    def plot_all(self):
        self.Filenames()
        self.load_all()
        sc.pl.spatial(self.adatas[1], color = "pct_counts_mt")
        # for i in range(len(self.adatas)):
        #     sc.pl.spatial(self.adatas[i], color = "pct_counts_mt", save = self.result_dir+f"\\IMG\\concatenation_GEO\\mtdata\\{i}.png", return_fig=False)



analise = Analysis()
analise.plot_tx()
