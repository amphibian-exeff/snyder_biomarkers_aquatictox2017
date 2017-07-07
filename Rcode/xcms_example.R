source("http://bioconductor.org/biocLite.R")
biocLite("xcms", dep=T)
biocLite("CAMERA")
library(xcms)
#change working directory to directory with .cdf files
setwd("/Library/Frameworks/R.framework/Versions/3.0/Resources/library/faahKO/cdf")
xset<-xcmsSet(method='centWave') #import cdf files
xset #check 
#retentio time correction
#xset1<-retcor(xset, method='obiwarp') 

xset2<-group(xset, bw=5, minfrac=0.5, mzwid=0.025)
# retention time correction and plots changes
xset3<-retcor(xset2, family="symmetric",plottype="mdevden" ) 

xset4<-group(xset3, bw=5, minfrac=0.5, mzwid=0.025)
xset5<-retcor(xset4, family="symmetric",plottype="mdevden" ) 

#fill peaks that are not really missing
xset3<-fillPeaks(xset4) #use xset that has had groups created

#create report of most significant differences
reporttab<-diffreport(xset3, "WT", "KO", "example", 10, metlin=0.15, h=480, w=640)

dr<-diffreport(xset3, filebase='variationA_vs_controlA', eicmax=100)
plotrt(xset1) #plots retention time plot
plotTIC(dr)



retcor(xset)
# example w faahKO data
library(multtest)
library(xcms)
library(faahKO)

setwd("/Library/Frameworks/R.framework/Versions/3.0/Resources/library/faahKO/cdf")

#plotting EIC example 
library(xcms)
library(faahKO)
xset<-xcmsSet(method='centWave') #import cdf files

xs<- group(xset)
xs<- fillPeaks(xs)
gt<-groups(xs)
eics <- getEIC(xs, mzrange=gt, rtrange = 200, groupidx = 1:nrow(gt))

pdf(file.path("tmp/%003d.pdf"), onefile = FALSE)
## you can also use png(file.path("tmp/%003d.png"), h=768, w=1024)
plot(eics, xs)
dev.off()

#plot individual EIC w example data
cdfpath <- ("~/Dropbox/amphib_dermalexposure/Biomarkers/Metabolomics/R_code/cdf")
cdffiles <- list.files(cdfpath, recursive = TRUE,full=T)
xr <- xcmsRaw(cdffiles[2])
xr
plotEIC(xr, m=c(550,551), rt=c(3200, 3800))
#plot individual EIC w tadpole data
cdfpath<-("~/Dropbox/amphib_dermalexposure/Biomarkers/Metabolomics/Marcia_140410_tadpoles_atrazine")
cdffiles <- list.files(cdfpath, recursive = TRUE,full=T)
xr <- xcmsRaw(cdffiles[4])
xr
plotEIC(xr,m=c(60,600), rt=c(421, 2900)) #rt = time range x-axis
plotEIC(xr,m=c(315,316), rt=c(421, 2900)) #rt = time range x-axis
# m is the ion 

#plot spectra
plotSpec(xr, ident=FALSE, vline=numeric(73))
tadpole1<-xcmsSet(cdffiles, method="centWave")
tadpole1_peaks<-findPeaks.centWave(xr)
plotPeaks(tadpole1_peaks)


