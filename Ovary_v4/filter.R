###!/usr/bin/env Rscript

# Suprime mensagens de inicialização dos pacotes
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

# Define uma função para identificar outliers
is_outlier <- function(x) {
  return(x > quantile(x, 0.75) + 1.5 * IQR(x))
}

# Define uma semente para garantir a reprodutibilidade
set.seed(1024)

# Define os IDs dos pacientes e os diretórios onde os dados estão armazenados
patient_id = c(1:8)
dir = list()

for (i in 1:length(patient_id)) {
  dir[[i]] = paste0(
    "~/sn_ovary/data/00_raw_data/output/SN-",
    patient_id[[i]] ,
    "/outs/per_sample_outs/SN-",
    patient_id[[i]],
    "/count/sample_filtered_feature_bc_matrix/"
  )
}

# Lê os dados de cada diretório e cria objetos Seurat, armazenando-os em uma lista
seurat.list <- lapply(dir, function(dir){
  mat <- Read10X(dir)
  seurat <- CreateSeuratObject(counts = mat, project = "sn_ovary", min.features = 100)
  return(seurat)
})

# Define os níveis dos fatores para cada paciente
levels = list()
for (i in 1:length(patient_id)){
 levels[[i]] = paste0("SN", patient_id[[i]])
}
names(seurat.list) = levels

# Adiciona a informação do paciente aos objetos Seurat
for (i in 1:length(patient_id)){
  seurat.list[[i]]$patient = paste0("SN", patient_id[[i]])
}

# Calcula a porcentagem de genes mitocondriais e identifica outliers para contagens e características de RNA
for (x in names(seurat.list)) {
  print(x)
  seurat.list[[x]]$percent.mt <-
    PercentageFeatureSet(seurat.list[[x]], pattern = "^MT-")
  seurat.list[[x]]@meta.data <-
    seurat.list[[x]]@meta.data %>% mutate(out_count = ifelse(is_outlier(nCount_RNA), "TRUE", "FALSE"))
  seurat.list[[x]]@meta.data <-
    seurat.list[[x]]@meta.data %>% mutate(out_feature = ifelse(is_outlier(nFeature_RNA), "TRUE", "FALSE"))
}

# Filtra os dados com base na porcentagem de genes mitocondriais, contagens de RNA e características de RNA
seurat.list.filt = list()
for (x in names(seurat.list)) {
  print(x)
  seurat.list.filt[[x]] <-
    subset(seurat.list[[x]],
           subset = percent.mt < 20 & out_count == FALSE &
             out_feature == FALSE)
}

# Conta o número de células antes e depois da filtragem
sum(unlist(lapply(seurat.list, ncol)))
sum(unlist(lapply(seurat.list.filt, ncol)))

# Define o diretório para salvar os objetos Seurat filtrados
path = "~/sn_ovary/data/01_preprocess/objects/"
for (x in 1:length(seurat.list.filt)) {
  print(x)
  saveRDS(object = seurat.list.filt[[x]],
          file = paste0(path,
                        names(seurat.list.filt[x]),
                        '.RDS',
                        collapse = NULL))
}

# Realiza a junção dos dados antes da filtragem
SN1 = seurat.list$SN1
seurat.list$SN1 = NULL
merge_prefilter = merge(SN1, seurat.list)

# Realiza a junção dos dados após a filtragem
SN1_f = seurat.list.filt$SN1
seurat.list.filt$SN1 = NULL
merge_posfilter = merge(SN1_f, seurat.list.filt)

# Produz gráficos de violino para visualizar a distribuição dos dados
p1 = VlnPlot(merge_prefilter, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), group.by = 'patient')
p2 = VlnPlot(merge_posfilter, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), group.by = 'patient')

# Salva os gráficos em um arquivo PNG
png("~/sn_ovary/data/01_preprocess/figures/Vln_QC.png", width = 14, height = 4, res = 300, units = 'in')
(p1 | p2) & xlab("")
dev.off()

# Produz gráficos de violino sem pontos para visualizar a distribuição dos dados
p1a = VlnPlot(merge_prefilter, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), group.by = 'patient', pt.size = 0)
p2a = VlnPlot(merge_posfilter, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), group.by = 'patient', pt.size = 0)

# Salva os gráficos sem pontos em um arquivo PNG
png("~/sn_ovary/data/01_preprocess/figures/Vln_QC_sempt.png", width = 14, height = 4, res = 300, units = 'in')
(p1a | p2a) & xlab("")
dev.off()
