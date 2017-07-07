#install.packages("MSeasy")
library(MSeasy) #used for clustering chromatograms
install.packages("Metab")
install.packages("xcms")

source("http://bioconductor.org/biocLite.R")
biocLite("mzR")
biocLite(c(“xcms”,”tcltk”,”R.utils”))
#source("http://bioconductor.org/biocLite.R")
#biocLite("xcms")
library(xcms)
# xcms is installed
#install package with demo data set 
install.packages("~/Downloads/faahKO_0.5.tgz", repos = NULL)

#### download metaXCMS http://metlin.scripps.edu/metaxcms/download.php

#### MALDIquant installed
library(MALDIquant)
demo(MALDIquant)
#MALDIquant works

#### for metadar installation https://code.google.com/p/metadar/wiki/Installation
install.packages(c("randomForest", "ROCR", "gplots", "beanplot", "subselect",
                   "glmnet", "mclust", "pROC", "plyr", "MASS"), dep=TRUE)

source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase", "simpleaffy"), dep=TRUE)

install.packages("gplots", dep=TRUE)
#download ihm.tar.gz and install
#download metadar.tar.gz and install
library(metadar)
# metadar works



