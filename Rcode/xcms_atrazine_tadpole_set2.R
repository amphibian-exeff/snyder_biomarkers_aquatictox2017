library(multtest)
library(xcms)
library(faahKO)

#### comparing control_methanol with highest dose
cdf.path<-("~/Dropbox/amphib_dermalexposure/Biomarkers/Metabolomics/Marcia_140410_tadpoles_atrazine/xcms_test")
list.files(cdf.path)

cdffiles <- list.files(cdf.path, recursive = TRUE, full.names = TRUE)

set2 <- xcmsSet(cdffiles)
#xset <- faahko
set2

set2 <- group(set2) #match peaks across samples

#correct for retention time offsets
set3 <- retcor(set2, family = "symmetric", plottype = "mdevden")

# group peaks again after correcting for time offset
set3 <- group(set3, bw = 10) 

#fill peaks
set4 <- fillPeaks(set3) #?


#normalize by peaks
set4_values<-groupval(set4)

#diff report
out.diff1 <- diffreport(set3, "c_meth", "d_1250", "example", 10, metlin = 0.15, h=480, w=640)

out.diff1 <- diffreport(set4_values, "c_meth", "d_1250", "example", 10, metlin = 0.15, h=480, w=640)
out.diff1[1:4,]

?xcmsRaw
out2<-xcmsRaw(cdffiles[1])
plotSpec(out2, vline=numeric(1790), rtrange=numeric(1700))

###########################################
###### comparing all doses
cdf.path<-("~/Dropbox/amphib_dermalexposure/Biomarkers/Metabolomics/Marcia_140410_tadpoles_atrazine/all_dose_xcms")
list.files(cdf.path)

cdffiles <- list.files(cdf.path, recursive = TRUE, full.names = TRUE)

set2 <- xcmsSet(cdffiles)
#xset <- faahko
set2

set2 <- group(set2) #match peaks across samples

#correct for retention time offsets
set3 <- retcor(set2, family = "symmetric", plottype = "mdevden")

# group peaks again after correcting for time offset
set3 <- group(set3, bw = 10) 

#fill peaks
set4 <- fillPeaks(set3) #?

#normalize by peaks
set4_values<-groupval(set4)

#diff report
out.diff1 <- diffreport(set4, "c_meth", "d_250", "250_meth", 10, metlin = 0.15, h=480, w=640)
out.diff2 <- diffreport(set4, "c_meth", "d_50", "50_meth", 10, metlin = 0.15, h=480, w=640)
out.diff2 <- diffreport(set4, "c_meth", "d_10", "10_meth", 10, metlin = 0.15, h=480, w=640)
out.diff2 <- diffreport(set4, "d_1250", "d_10", "d10_1250", 10, metlin = 0.15, h=480, w=640)
out.diff2 <- diffreport(set4, "c_meth", "d_1250", "1250_meth", 10, metlin = 0.15, h=480, w=640)

#out.diff1 <- diffreport(set4_values, "c_meth", "d_1250", "example", 10, metlin = 0.15, h=480, w=640)
out.diff1[1:4,]





