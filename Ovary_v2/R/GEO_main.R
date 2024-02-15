########
#Stage = 3C
#Asian
#Age = 49
#Response to nact = Excellent
#GEO
#########
library(STdeconvolve)
library(SpotClean)
movary_raw = read10xRaw("data/GEO1/raw_feature_bc_matrix/")
dim(movary_raw)
head(movary_raw)

slide_info = read10xSlide("data/GEO1/spatial/tissue_positions_list.csv",
                          "data/GEO1/spatial/II21472_tissue_lowres_image.png",
                          "data/GEO1/spatial/scalefactors_json.json")
# slide
slide = createSlide(movary_raw,
                    slide_info)
visualizeSlide(slide)


#criando countmat na mao
library(Matrix)
countmat = Matrix::readMM("data/GEO1/raw_feature_bc_matrix/matrix.mtx.gz")
dim(countmat)
head(countmat)

barcode = read.csv("data/GEO1/raw_feature_bc_matrix/barcodes.tsv.gz", header = FALSE)
dim(barcode)

gene = read.csv("data/GEO1/raw_feature_bc_matrix/features.tsv.gz", header = FALSE, sep="\t")
dim(gene)
head(gene)
View(gene)

# COUNTMAT
colnames(countmat) <- barcode[,1]
rownames(countmat) <- make.unique(gene[,2])

head(countmat)

# POS
pos.info = read.csv("data/GEO1/spatial/tissue_positions_list.csv", header = FALSE)
head(pos.info)
pos = pos.info[,c("V6", "V5")]
dim(pos)
row.names(pos) = pos.info[,1]
head(pos)
plot(pos, pch = 16)

# fix order position
dim(pos)
pos = pos[colnames(countmat),]
pos[,2] = -pos[,2]
pos[,2] = pos[,2] -min(pos[,2])
head(pos)
dim(pos)
dim(countmat)
plot(pos)


# position done

plot(pos, pch = 16)

### genxp
genxp = colSums(countmat)
head(genxp)
g = c("APOE", "LYZ")
g = "FOLR2"
g = "LYZ"
g = "TMEM125"
gen = countmat[g,]
par(mfrom= c(1,1))
MERINGUE::plotEmbedding(pos, col = genxp)
MERINGUE::plotEmbedding(pos, col = gen)



# filtros

corpus <- restrictCorpus(countmat, plot = TRUE)
dim(corpus)

library(Matrix)
## choose optimal number of cell-types
ks <- seq(from = 2, to = 18, by = 1) # range of K's to fit LDA models with given the input corpus
ldas <- fitLDA(as.matrix((corpus),
                         Ks = ks,
                         ncores = parallel::detectCores(logical = TRUE) - 1, # number of cores to fit LDA models in parallel
                         plot=TRUE, verbose=TRUE))


## get best model results
optLDA <- optimalModel(models = ldas, opt = "8")


## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0)
deconProp <- results$theta
deconGexp <- results$beta

head(deconGexp)
head(deconProp)
dim(deconProp)
dim(pos)
dim(pos[rownames(deconProp),])
## visualize deconvolved cell-type proportions
head(rownames(deconProp))
head(pos)
colnames(pos) <- c('x', 'y')
head(pos)

#error
pos <- pos[rownames(deconProp)]
head(ppos)
vizAllTopics(deconProp, pos, r = 2.5)





