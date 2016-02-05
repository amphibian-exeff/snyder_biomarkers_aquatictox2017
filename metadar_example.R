## load metadar library
library("metadar")
#### there are different ways to read the data
## Data consists of three components
## (1) data
## (2) sample meta data
## (3) variable meta data
#### Data is represented as a Dataset object

### one input file (CSV file) example
dat1 <- new("Dataset", "examples/egData_standalone.csv")
### two input files (CSV files)
# example 1
dat <- new("Dataset", "examples/egData1.csv",
           "examples/egSampleMetadata.csv")
# example 2
dat2 <- new("Dataset", "examples/egData2.csv",
            "examples/egSampleMetadata.csv")

### three files (CSV file) example
dat3 <- new("Dataset", "examples/egData1.csv",
            "examples/egSampleMetadata.csv",
            "examples/egVarMetaData.csv")

### Filter the data by removing compounds that have too many zero's.
### By minimum number of non-zero values
dat <- zeroFiltering(dat, minNfound=5)
### By minumim number of non-zero values per group defined by the covariate
dat <- zeroFiltering(dat, minNfound=5, covariate="Class")
### By minimum percentage of non-zero values
dat <- zeroFiltering(dat, pctNfound=70)
### By minumim percentage of non-zero values per group defined by the covariate
dat <- zeroFiltering(dat, pctNfound=70, covariate="Phenotype")

### zero imputation
dat <- zeroImputation(dat)
# OR
dat <- zeroImputation(dat, covariate="DiagnosticGroup")

### log transform the data
dat <- log2(dat)
# OR
dat <- log10(dat)
# OR
dat <- log(dat, base=2)

#### Here is an example of some simple data analysis of
#### one data set which had 24 samples from
#### two phenotypes measured under two time points
#### we are interested in the difference between
#### the case versus control phenotype at each time point
#### I am assuming that the dataset is called 'dat'
#### and that there are columns titled 'Phenotype'
#### and 'Time' in sample meta data

### take a subset of the data (in this case we want to analyze each time point separately)
dat.3d <- get.array.subset(dat, "Time", "3d")
dat.24h <- get.array.subset(dat, "Time", "24h")


### Perform t-test on 3d data
tt.3d <- univariateTTest(dat.3d, "Phenotype")

### Paired t-test on 3d data
tt.3d <- univariateTTest(dat.3d, "Phenotype", paired=TRUE)

### Fold change on 3d data
fc.3d <- foldChange(dat.3d, "Phenotype", is.logged=FALSE)
# OR
fc.3d <- foldChange(dat.3d, "Phenotype", is.logged=TRUE, log.base=2)

### Fold change when samples are paired!
fc.3d <- foldChange(dat.3d, "Phenotype", paired=TRUE)
# OR
fc.3d <- foldChange(dat.3d, "Phenotype", paired=TRUE, is.logged=TRUE, log.base=2)

### Compute correlation coefficient of each variable versus one meta data variable
corrs <- univariateCorrelation(dat, covariate="Fins")
# OR
corrs <- univariateCorrelation(dat, covariate=fins) ## where fins is the numerical variable
# OR
corrs <- univariateCorrelation(dat, covariate="Fins", method="pearson")
# OR
corrs <- univariateCorrelation(dat, covariate=fins, method="pearson")

### odds ratio of each variable with respect to a binary outcome
aucs <- univariateOddsRatio(dat, "Phenotype")

### The area under the ROC curve (AUC) of each variable
aucs <- univariateAUC(dat, "Phenotype")

### Hierarchical clustering
plot(colhclust(dat, labels="SampleName", color="Phenotype"))

### one way anova test
aovres <- oneWayAnova(dataset, covariate="Class")

#### generate a heatmap from one way anova result
log2foldchanges <- log2(aovres[, grep("ratio", colnames(aovres))])
rownames(log2foldchanges) <- rownames(aovres)
pvalues <- aovres[, grep("p adj", colnames(aovres))]
rownames(pvalues) <- rownames(aovres)
take <- apply(pvalues, function(x) any(x < 0.05))

ihm(log2foldchanges[take,], pvalues[take,], clusterColumns=TRUE, clusterRows=TRUE,
    device="pdf", file="Heatmap_sig.pdf", cellnotesizeCorrection=2)

### two way anova
aov2 <- twoWayAnova(dat, covariate1="Drug", covariate2="Time")

### draw a heatmap with the two way anova result
ihm(aov2[["RatiosTable"]][, grep("ratio", aov2[["RatiosTable"]])], aov2[["PvalsTable"]])

### model based clustering is explained at
### URL: http://code.google.com/p/metadar/wiki/Mclustering

## PCA
pc <- new("PCA", dataset, center=TRUE, scale=FALSE)
pc$scoresplot(annotation="Class", color="Class", pcs=1:2, cex=2, labelpos=2)
pc$loadingsplot(pcs=1:2, labelpos=4, cex=2)

