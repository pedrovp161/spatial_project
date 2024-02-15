####
#
#SEURAT
#
####
library(SeuratData)


barcode = read.csv("data/raw_feature_bc_matrix/barcodes.tsv.gz", header = FALSE)
View(barcode)
dim(barcode)

countmat = Matrix::readMM("data/raw_feature_bc_matrix/matrix.mtx.gz")
head(countmat)
dim(countmat)
# COUNTMAT
colnames(countmat) <- barcode[,1]
rownames(countmat) <- make.unique(gene[,2])
head(countmat)
dim(countmat)



Slice1 = Read10X_Image("data/spatial/GSM5708485_II21472_tissue_lowres_image.png.gz")

gene = read.csv("data/raw_feature_bc_matrix/features.tsv.gz", header = FALSE, sep="\t")
head(gene)
dim(gene)


### datamodel
brain <- LoadData("stxBrain", type = "anterior1")






library(Seurat)

obj = CreateSeuratObject(counts = countmat,
                         assay = "Spatial",
                         slice = "slice1",
                         filter.matrix = TRUE,
                         to.upper = FALSE,
                         image = Slice1)



