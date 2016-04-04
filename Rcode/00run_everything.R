library(reshape2) # for reformatting table

gt_root <- "~/git/Snyderetal2016_biomarkers/"
gt_data <- paste(gt_root, "data_in/", sep="")
gt_dir <- paste(gt_data, "gt_stage22_TOF041515/xcms_polar/results/", sep="")


source(gt_root, "01boxplots_pathways.R", sep="")