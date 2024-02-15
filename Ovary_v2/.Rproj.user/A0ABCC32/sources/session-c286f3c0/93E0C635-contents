
library(Seurat)

Idents(ov_combined) <- ov_combined$seurat_clusters

pdf("./figures/21480_HE_Cluster_060221.pdf")
SpatialDimPlot(ov_combined, images = "PR_2", pt.size.factor = 3) + 
    theme(legend.position = "right")
dev.off()


pdf("./figures/21481_HE_Cluster_060221.pdf")
SpatialDimPlot(ov_combined, images = "PR_3", pt.size.factor = 3) + 
    theme(legend.position = "right")
dev.off()
