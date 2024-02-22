from utils.leitura import Leitura
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
import os

# definindo seed para reprodutibilidade
random.seed(10)
print(f"seed utilizada: {random.random()}")

# Define função para identificação de outliers.
def is_outlier(x):
    return x > np.quantile(x, 0.75) + 1.5 * (np.percentile(x, 75) - np.percentile(x, 25))

# define função para contagem de células em todas as amostras
def check_cell_sum(diretorio):
    somatorio_de_celulas = []
    for i in diretorio:
        adata = diretorio[i]
        cells_before_by_sample = sum([adata.n_obs for adata in diretorio[i]])
        somatorio_de_celulas.append(cells_before_by_sample)
    return print(sum(somatorio_de_celulas)) # falta o antes

# definindo diretorio
DIR = os.path.dirname(__file__)

# Objetificando a classe e passando diretório
ler = Leitura(DIR)

# Obtendo diretório com dados
adatas_dir = ler.listar_raw_GEO_dic_com_dados()

# Calcula a porcentagem de genes mitocondriais e identifica outliers para contagens e características de RNA

for i in adatas_dir:
    adata = adatas_dir[i]
    print(adata)
    adata.var["mt"] = adata.var["gene_ids"].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata.obs['pct_counts_mt'] = np.sum(adata[:, adata.var_names.str.startswith('MT-')].X, axis=1) / np.sum(adata.X, axis=1) * 100 # pct_counts_mt é a porcentagem de reads mitocondriais e é equivalente 
    adata.obs['out_count'] = is_outlier(adata.obs['total_counts']).astype(str) # total_counts é o numero de reads e é o equivalente a n_counts
    adata.obs['out_feature'] = is_outlier(adata.obs['n_genes_by_counts']).astype(str) # n_genes_by_counts = n_genes

    # Filtra os dados com base na porcentagem de genes mitocondriais, contagens de RNA e características de RNA.
    adata = adata[adata.obs["pct_counts_mt"] < 20]
    print(adata)


### Contagem do número de células antes e depois da filtragem
# check_cell_sum(diretorio=adatas_dir)
# check_cell_sum(diretorio=adatas_dir_raw)

# Define o diretório para salvar os objetos Anndata filtrados
