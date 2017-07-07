library(xcms); library(faahKO); data(faahko)

tapply(peaks(faahko)[,"into"], as.factor(peaks(faahko)[,"sample"]), sum)

tapply(peaks(faahko)[,"into"], as.factor(peaks(faahko)[,"sample"]),
       sum)

peakTable <- function(xs){
  if (nrow(xs@groups) > 0) {
    groupmat <- groups(xs)
    ts <- data.frame(cbind(groupmat,groupval(xs, "medret",
                                             "into")),row.names = NULL)
    cnames <- colnames(ts)
    colnames(ts) <- cnames
  } else if (length(xs@sampnames) == 1)
    ts <- xs@peaks
  else stop ('First argument must be a xcmsSet with group information
or contain only one sample.')
  ts
}


cdf.path<-setwd("~/Dropbox/amphib_metabolomics/DATA/toad_stage22_round1/polar/polar_toad1_cdf/xcms")
list.files(cdf.path)

cdffiles <- list.files(cdf.path, recursive = TRUE, full.names = TRUE)

set2 <- xcmsSet(cdffiles)
xs <- xcmsSet(set2)
xg <- group(xs)
pt <- peakTable(xg)
write.table(pt,file="peaktable.csv",row.names=F)
######
xcmsSet.TIC<-function(object,...){
  if(class(object) == "xsAnnotate"){
    object<-object@xcmsSet
  }
  Intensity<-rowMeans(groupval(object, "medret", "maxo"))
  rt<-object@groups[,"rtmed"]
  idx<-order(rt)
  
  plot(x=rt[idx], y=Intensity[idx], type="b",
       ylab="Intensity", xlab="Retention time (sec)",
       main="Peak picked TIC", ...)
}
xcmsSet.TIC(xg)