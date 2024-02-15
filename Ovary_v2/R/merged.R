library(Matrix)
list = c("data/GEO1/raw_feature_bc_matrix/matrix.mtx.gz", "data/GEO2/raw_feature_bc_matrix/matrix.mtx.gz",
         "data/GEO3/raw_feature_bc_matrix/matrix.mtx.gz", "data/GEO4/raw_feature_bc_matrix/matrix.mtx.gz",
         "data/GEO5/raw_feature_bc_matrix/matrix.mtx.gz", "data/GEO6/raw_feature_bc_matrix/matrix.mtx.gz")
matrix_list = list()
b = 0
for (i in list){
  b = b+1
  matrix_info = Matrix::readMM(i)
  matrix_list[b] = matrix_info
}

list = c("data/GEO1/raw_feature_bc_matrix/barcodes.tsv.gz", "data/GEO2/raw_feature_bc_matrix/barcodes.tsv.gz",
         "data/GEO3/raw_feature_bc_matrix/barcodes.tsv.gz", "data/GEO4/raw_feature_bc_matrix/barcodes.tsv.gz",
         "data/GEO5/raw_feature_bc_matrix/barcodes.tsv.gz", "data/GEO6/raw_feature_bc_matrix/barcodes.tsv.gz")
barcode_list = list()
b=0
for (i in list) {
  b=b+1
  barcode = read.csv(i, header = FALSE)
  barcode_list[b] = barcode
}


# COUNTMAT
gene = read.csv("data/GEO1/raw_feature_bc_matrix/features.tsv.gz", header = FALSE, sep="\t")
for (i in 1:6) {
  colnames(matrix_list[[i]]) <- barcode_list[[i]]
  rownames(matrix_list[[i]]) <- make.unique(gene[,2])
  countmat_E = matrix_list
}




#pos
list = c("data/GEO1/spatial/tissue_positions_list.csv", "data/GEO2/spatial/tissue_positions_list.csv",
         "data/GEO3/spatial/tissue_positions_list.csv", "data/GEO4/spatial/tissue_positions_list.csv",
         "data/GEO5/spatial/tissue_positions_list.csv", "data/GEO6/spatial/tissue_positions_list.csv")
pos_list = list()
pos_barcodes = list()
for (i in list) {
  pos.info = read.csv(i, header = FALSE)
  pos_list = append(pos_list, pos.info[,c(5,6)], after = 2)
  pos_barcodes = append(pos_barcodes, as.data.frame(pos.info[,1]))
}

# pos dos Exelent_responders
pos_1 = cbind(as.data.frame(pos_list[[1]]), as.data.frame(pos_list[[2]]))
pos_2 = cbind(as.data.frame(pos_list[[3]]), as.data.frame(pos_list[[4]]))
pos_3 = cbind(as.data.frame(pos_list[[5]]), as.data.frame(pos_list[[6]]))
pos_4 = cbind(as.data.frame(pos_list[[7]]), as.data.frame(pos_list[[8]]))
pos_5 = cbind(as.data.frame(pos_list[[9]]), as.data.frame(pos_list[[10]]))
pos_6 = cbind(as.data.frame(pos_list[[11]]), as.data.frame(pos_list[[12]]))

pos_E = list(pos_1,pos_2,pos_3,pos_4,pos_5, pos_6)
x=1
barcode = as.data.frame(pos_barcodes[[2]])
while (x<=length(pos_E)) {
  colnames(pos_E[[x]]) = c("x","y")
  row.names(pos_E[[x]]) = barcode[,1] # add barcodes
  pos_E[[x]] = pos_E[[x]][colnames(countmat_E[[x]]),]
  pos_E[[x]] = rev(pos_E[[x]])
  pos_E[[x]][,2] = -pos_E[[x]][,2]
  pos_E[[x]][,2] = pos_E[[x]][,2] -min(pos_E[[x]][,2])
  x=x+1
}

# plots
library(MERINGUE)

MERINGUE::plotEmbedding(pos_E[[5]], col = colSums(countmat_E[[5]]))

g = "FOLR2"
gen = countmat_E[[5]][g,]
MERINGUE::plotEmbedding(pos_E[[5]], col = gen)

#results
save(pos_E, countmat_E, file = "results/mergeddata1")


















#countmatmerge = merge(as.matrix(countmat_list[[1]]), as.matrix(countmat_list[[2]]))






# library(SeuratObject)
# supermat = RowMergeSparseMatrices(countmat_list[[1]], countmat_list[[2]]) 






