import scanpy as sc
import anndata
import os
print(os.getcwd())
sc.set_figure_params(figsize=(6,6), frameon=False)

adata = sc.read_10x_h5(f"{os.getcwd()}/data/5k_human_ovarian_tumor_CNIK_5pv2_filtered_feature_bc_matrix.h5")
adata.obs["sample"] = "10x1"

print(adata.var_names)