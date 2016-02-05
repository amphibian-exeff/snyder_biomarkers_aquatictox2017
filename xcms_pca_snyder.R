source("http://bioconductor.org/biocLite.R")
biocLite("xcms")
biocLite("pcaMethods")
library(pcaMethods)
library(xcms)

if(Sys.info()[4]=="DZ2626UTPURUCKE"){
  cdf.path <- "k:\\git\\Snyderetal2016_biomarkers\\cdf"
}
if(Sys.info()[4]=="marcias laptop name here"){
  cdf.path<-("~/Dropbox/amphib_dermalexposure/Biomarkers/Metabolomics/R_code/cdf")
}

list.files(cdf.path)
cdffiles <- list.files(cdf.path, recursive = TRUE,full=T)

#http://metabolomics-forum.com/viewtopic.php?f=26&t=124

#construction of xcmsSet objects. It finds peaks in batch mode and 
#pre-sorts files from subdirectories into different classes suitable for grouping
faahko <- xcmsSet(cdffiles)
#grouping (or alignment) methods
#density-based approach described by Colin Smith (Scripps) et al (2006) 
#returns an xcmsSet object with peak group assignments and statistics.
faahko<-group(faahko, method = "density")
#For each sample, identify peak groups where that sample is not represented. 
#For each of those peak
#groups, integrate the signal in the region of that peak group and create a new peak.
faahko<-fillPeaks(faahko)
#Generate a matrix of peak values with rows for every group and columns for 
#every sample. The value included in the matrix can be any of the columns from 
#the xcmsSet peaks slot matrix. Collisions
#where more than one peak from a single sample are in the same group get resolved
values<-groupval(faahko, intensity="into")

## Get peak intensity values
#??values <- groupval(values, intensity="into")

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
plotPcs(pca.result, type = "scores", col=as.integer(sampclass(faahko)) + 1)
dev.new()
## look to see if the model is valid

plot(pca.result)

#Then have a look with the validity of the model using a Q2 values.
#Loadings will let you know which metabolite is effecting the model 
#and which metabolites are causing the differentiation.
plot(loadings[,1], loadings[,2], pch=16, cex=0.5)
text(loadings[,1], loadings[,2], rownames(values), col="red", cex=0.5)

## MDS Example
library(MASS)

## MDS, For the MDS example we'll use the same data as before.
data.dist <- dist(data)
mds <- isoMDS(data.dist)
plot(mds$points, type = "n")
text(mds$points, labels = rownames(data),col=as.integer(sampclass(faahko)) + 1)
