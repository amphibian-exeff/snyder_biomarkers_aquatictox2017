library(xcms)
source("http://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
library(pcaMethods)

cdf.path<-("~/Dropbox/amphib_dermalexposure/Biomarkers/Metabolomics/R_code/cdf")
list.files(cdf.path)
cdffiles <- list.files(cdf.path, recursive = TRUE,full=T)
faahko <- xcmsSet(cdffiles)
faahko<-group(faahko)
faahko<-fillPeaks(faahko)
values<-groupval(faahko, intensity="into")

## Get peak intensity values
values <- groupval(values, intensity="into")

## Numerical matrix with with samples in rows and variables as columns
data <- t(values)

## Normalize each mass signal to max=1
for (r in 1:ncol(data))
  data[,r] <- data[,r] / max(data[,r])

##  PCA Example
pca.result <- pca(data, nPcs = 3, scale="none", cv="q2")

## Get the estimated principal axes (loadings)
loadings <- pca.result@loadings

## Get the estimated scores
scores <- pca.result@scores

## Now plot the scores
plotPcs(pca.result, type = "scores", col=as.integer(sampclass(xsg)) + 1)
dev.new()
## look to see if the model is valid
plot(pca.result)