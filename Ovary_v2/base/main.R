genxp = t(countmat)
head(genxp)

colSums(genxp)

hist(log10(colSums(genxp>0)+1))
nrow(genxp)

"MS4A1" %in% colnames(genxp)

g = "MS4A1"

