library(xcms)
library(faahKO)
library(pls)

cdf.path<-("~/Dropbox/amphib_dermalexposure/Biomarkers/Metabolomics/R_code/cdf")
list.files(cdf.path)
cdffiles <- list.files(cdf.path, recursive = TRUE,full=T)
faahko <- xcmsSet(cdffiles)
faahko<-group(faahko)
faahko<-fillPeaks(faahko)
values<-groupval(faahko)

val<-data.frame(as.numeric(sampclass(faahko)))
val[,2]<-round(t(values),4)
colnames(val)<-c("class", "Met")

faahko.pls <- plsr(class ~ Met, ncomp = 3, data=val, validation = "LOO")
summary(faahko.pls)
plot(RMSEP(faahko.pls), legendpos = "topright") #Root Mean Square Error of Prediction = RMSEP
dev.new()
plot(faahko.pls, ncomp = 2, asp = 1, line = TRUE)
dev.new()
## to do this manually

scoreplot(faahko.pls, comps = 1:3, identify = FALSE, type = "p", col=rep(c("red", "blue"), c(6,6)), pch=16 )
loadingplot(faahko.pls, comps = 1:3, identify = TRUE, type= "p")
## use identify == TRUE if you want to identify your loadings
## these will be printed to the R console

summary(faahko.pls)
plot(faahko.pls, ncomp = 3)
plot(faahko.pls, "scores", type="p", col=rep(c("red", "blue")))
plot(faahko.pls, "loadings", comps = 1:3)
plot(faahko.pls, "coef", ncomp = 3)
plot(faahko.pls, "val")
plot(faahko.pls, "val", val.type = "MSEP", estimate = "CV")
