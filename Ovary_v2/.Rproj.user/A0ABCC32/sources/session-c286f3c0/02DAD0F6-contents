library(BiocManager)
library(SpotClean)
library(STdeconvolve)

movary_raw <- read10xRaw("data/10x1/raw_feature_bc_matrix/")

dim(movary_raw)
head(movary_raw)

library(Matrix)
countmat <- Matrix::readMM("data/10x1/raw_feature_bc_matrix/matrix.mtx.gz")
dim(countmat)
head(countmat)

barcode.info = read.csv("data/10x1/raw_feature_bc_matrix/barcodes.tsv.gz", header=FALSE)

barcode <- barcode.info[,1]


gene.info <- read.csv("data/10x1/raw_feature_bc_matrix/features.tsv.gz", header=FALSE, sep="\t")
gene.info


gene <- gene.info[,2]
colnames(countmat) <- barcode
rownames(countmat) <- make.unique(gene)

pos.info <- read.csv('data/10x1/spatial/tissue_positions_list.csv', header=FALSE)
pos <- pos.info[,c(5,6)]
rownames(pos) <- pos.info[,1]
head(pos)
plot(pos)

## fix order
pos <- pos[colnames(countmat),]




## fix position
pos <- pos[, c(2,1)]
pos[,2] <- -pos[,2] 
pos[,2] <- pos[,2] - min(pos[,2])
head(pos)

library(Matrix)
libsize <- colSums(countmat)
par(mfrow=c(1,1))
MERINGUE::plotEmbedding(pos, col=libsize)

g <- 'FOLR2'
g <- "Foxl2"
g <- 'LYZ'
gexp <- countmat[g,]
gexp
MERINGUE::plotEmbedding(pos, col=gexp[rownames(pos)], main=g)


dim(pos)
head(pos)
plot(pos)

#### run STdeconvolve

## remove pixels with too few genes
counts <- cleanCounts(countmat,
                      min.lib.size = 10^3.75,
                      min.reads = 10^3.5,
                      plot = TRUE)
dim(counts)
par(mfrow=c(1,1))
MERINGUE::plotEmbedding(pos[colnames(counts),], col=colSums(counts))

g %in% rownames(counts)

## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
dim(corpus)
head(corpus)
library(Matrix)
## choose optimal number of cell-types
ks <- seq(from = 2, to = 18, by = 1) # range of K's to fit LDA models with given the input corpus
ldas <- fitLDA(as.matrix((corpus),
                         Ks = ks,
                         ncores = parallel::detectCores(logical = TRUE), # number of cores to fit LDA models in parallel
                         plot=TRUE, verbose=TRUE))

#ldas <- fitLDA((as.matrix(corpus)), Ks =25)


## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")


## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
head(deconProp)
dim(deconProp)
dim(pos)
dim(pos[rownames(deconProp),])
## visualize deconvolved cell-type proportions
ppos <- pos[rownames(deconProp),]
colnames(ppos) <- c('x', 'y')
vizAllTopics(deconProp, ppos, r=2.5)
