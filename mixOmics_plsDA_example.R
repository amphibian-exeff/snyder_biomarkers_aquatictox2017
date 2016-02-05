library(mixOmics)

data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4]) # Y is a factor, we chose it as the
# time points of necropsy              

## PLS-DA function
result <- plsda(X, Y, ncomp = 3) # where ncomp is the number of components wanted

## sPLS-DA function
result <- splsda(X, Y, ncomp = 3, keepX = c(10, 10, 10)) # where keepX is the number
# of variables selected
# for each components