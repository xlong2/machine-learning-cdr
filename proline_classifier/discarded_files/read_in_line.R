# read in the spline of what the training set is
library("FactoMineR")
data("decathlon")
res.pca <- PCA(decathlon, quanti.sup = 11:12, quali.sup = 13)


