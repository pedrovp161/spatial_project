
library(tidyverse)
library(Seurat)
library(ggpubr)

VlnPlot(ov_combined, features = c("TP53","BRCA1","BRCA2","PTEN"), split.by = 'Group', ncol = 1, 
        slot = 'scale.data')

exps <- GetAssayData(ov_combined, assay = 'SCT', slot = 'data')

dim(exps)

cancer_exps <- cbind(ov_combined@meta.data, exps[c("TP53","BRCA1","BRCA2","PTEN"),] %>%
    as.matrix() %>%
    t() %>%
    as_tibble()) %>%
    as_tibble()


pdf('./figures/Cancer_Gene_Expression_042121.pdf', width = 15, height = 7)
ggviolin(cancer_exps, fill = 'Group', y = 'TP53', x = 'seurat_clusters', palette = 'Set2') +
    stat_compare_means(aes(group = Group), label = "p.format")
ggviolin(cancer_exps, fill = 'Group', y = 'BRCA1', x = 'seurat_clusters', palette = 'Set2') +
    stat_compare_means(aes(group = Group), label = "p.format")
ggviolin(cancer_exps, fill = 'Group', y = 'BRCA2', x = 'seurat_clusters', palette = 'Set2') +
    stat_compare_means(aes(group = Group), label = "p.format")
ggviolin(cancer_exps, fill = 'Group', y = 'PTEN', x = 'seurat_clusters', palette = 'Set2') +
    stat_compare_means(aes(group = Group), label = "p.format")
dev.off()
