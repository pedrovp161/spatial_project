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

def is_outlier_low(x):
    return x < np.quantile(x, 0.25) - 1.5 * (np.percentile(x, 75) - np.percentile(x, 25))

# define função para contagem de células em todas as amostras
def check_cell_sum(dicionario):
    somatorio_de_celulas = []
    for i in dicionario:
        adata = dicionario[i]
        cells_before_by_sample = sum([adata.n_obs for adata in dicionario[i]])
        somatorio_de_celulas.append(cells_before_by_sample)
    return print(sum(somatorio_de_celulas)) # falta o antes

# definindo função pra processamento
def processar(adatas_dir):
    for i in adatas_dir:
        adata = adatas_dir[i]
        adata.var["mt"] = adata.var["gene_ids"].str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        # marcação acima do 3 quartil
        adata.obs['out_counts'] = is_outlier(adata.obs['total_counts']).astype(str) # total_counts é o numero de reads e é o equivalente a n_counts

        #marcação e corte antes do 1 quartil
        adata.obs["out_counts_low"] = is_outlier_low(adata.obs["total_counts"]).astype(str)
        adata.obs["total_counts"] = adata.obs["total_counts"][~adata.obs["out_counts_low"].astype(bool)]

        # outliers de genes marcados
        adata.obs['out_genes'] = is_outlier(adata.obs['n_genes_by_counts']).astype(str) # n_genes_by_counts = n_genes

        # Filtra os dados com base na porcentagem de genes mitocondriais, contagens de RNA e características de RNA.
        adata = adata[adata.obs["pct_counts_mt"] < 20]
        adatas_dir[i] = adata

# definindo diretorio
DIR = os.path.dirname(__file__)

# Objetificando a classe e passando diretório
ler = Leitura(DIR)

# Obtendo diretório com dados
adatas_dir = ler.listar_raw_GEO_dic_com_dados()

# salvando dado raw
adatas_dir_raw = adatas_dir

# Calcula a porcentagem de genes mitocondriais e identifica outliers para contagens e características de RNA
processar(adatas_dir)


### Contagem do número de células antes e depois da filtragem
check_cell_sum(dicionario=adatas_dir)
check_cell_sum(dicionario=adatas_dir_raw)

# Define o diretório para salvar os objetos Anndata filtrados
output_dir = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__),"..")), "Ovary_v4", "data", "Resultados", "Filtro_geral")

# Verifica se o diretório de saída existe, se não, cria
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for nome, adata in adatas_dir.items():
    output_file_path = os.path.join(output_dir, f"{nome}.h5ad")
    sc.write(output_file_path, adata)





