library(reshape2) # for reformatting table

gt_root <- "~/git/Snyderetal2016_biomarkers/"
gt_rcode <- paste(gt_root, "Rcode/", sep = "")
gt_data <- paste(gt_root, "data_in/", sep= "")
gt_figures <- paste(gt_root, "figures/", sep="")
gt_dir <- paste(gt_data, "gt_stage22_TOF041515/xcms_polar/results/", sep="")


source(paste(gt_rcode, "01boxplots_pathways.R", sep=""))
