#### gather data GT ####
# import metabolite names for IDs 
#ID_METABOLITE_MattFinalID2_112315 <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/ID_METABOLITE_MattFinalID2_112315.csv")
#View(ID_METABOLITE_MattFinalID2_112315)
gt_metabolites <- paste(gt_dir, "ID_METABOLITE_GTMattFinalID2_120115.csv", sep = "")
file.exists(gt_metabolites)
ID_METABOLITE_GTMattFinalID2_120115 <- read.csv(gt_metabolites)
#View(ID_METABOLITE_GTMattFinalID2_120115)

# import pre-processed bin data 
file_data_class <- paste(gt_data, "class.csv", sep="")
file.exists(file_data_class)
class <- read.csv(file_data_class)
file_processed_table <- paste(gt_data, "ProcessedTable.csv", sep="")
file.exists(file_processed_table)
ProcessedTable <- read.csv(file_processed_table)
data2<-cbind(class$V1, ProcessedTable)
#View(data2)
#subset data by the 200 top ranked bins with SVM-RFE
data2.sub <- data2[,names(data2) %in% ID_METABOLITE_GTMattFinalID2_120115$peaks3 ]
dim(data2.sub)
# add class and name column back
head(names(data2))
data2.sub2<-cbind(data2$X,class$V1, data2.sub)
#View(data2.sub2)
# transform data2 from wide form to long form using lib(reshape2)
data.reshape2<-melt(data2.sub2, id.vars=c("class$V1", "data2$X"))
# use apply to add column where bin name = metabolite ID
dim(data.reshape2)
#add column for metabolite ID based on bin_name
data.reshape2$metabolite<-c("unknown")
dim(data.reshape2)

#View(data.reshape2)
names(data.reshape2)
x<-data.reshape2

class(data.reshape2)
#rename retention time mass fragment bins with metabolite ID with merge
x2<-merge(x,ID_METABOLITE_GTMattFinalID2_120115, by.x="variable", by.y="peaks3" )
#View(x2)
unique(x2$mattID1)
#drop bins with similiarity value <700
x3<-x2[x2$sim1>699,]
dim(x3)

#### boxplots GT ####
names(x3)
head(x3)
#boxplot(x3$value~x3$"class$V1")
# box plots for all doses TOF ID
unique.metabolites<-unique(x3$mattID1)
for (metabolite in unique.metabolites){
  temp.data<-x3[which(x3$mattID1==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", main=paste(metabolite))
}
#subset dataset to control and highest dose
xc1250<-subset(x3,x3$`class$V1`=="1" | x3$`class$V1`=="5" )
#View(xc1250)
dim(xc1250)
# boxplot for control and 1250 dose TOF ID
unique.metabolites<-unique(x3$mattID1)
for (metabolite in unique.metabolites){
  temp.data<-xc1250[which(xc1250$mattID1==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", main=paste(metabolite))
}

#rename metabolites based on KEGG names
#import csv for name translation
gt_id_kegg <- paste(gt_data, "GT_ID_Kegg_120215.csv", sep="")
GT_ID_Kegg_120215 <- read.csv(gt_id_kegg)
names(GT_ID_Kegg_120215)
x4<-merge(x3, GT_ID_Kegg_120215, by.x = "mattID1", by.y = "mattID1" )
#test<-x3[which(x3$mattID1=="n-acetyl-l-tyrosine"),] #ys
names(x4)
#View(x4)
dim(x4)
#subset data for control and 1250 dose w Kegg ID
xc1250b<-subset(x4,x4$`class$V1`=="1" | x4$`class$V1`=="5" )
# print pdf of all boxplots
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
gt_boxplots <- paste(gt_figures, "GT_boxplots120315.pdf", sep="")
pdf(gt_boxplots)
# boxplot for control and 1250 dose TOF ID
unique.metabolites<-unique(xc1250b$KeggMatch)
for (metabolite in unique.metabolites){
  temp.data<-xc1250b[which(xc1250b$KeggMatch==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", main=paste(metabolite))
}
dev.off()
##################################################################################
#### plot individual box plots for metabolites in top pathways from overlap GT and AT
# box plots for aminoacyl-trna biosynthesis GT
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots")
setwd(gt_figures)
#pdf("GT_boxplots_aminoacyl_120315.pdf")
# in aminoacyl tRNA biosynthesis pathway
aminoacyl.metabs.gt<-as.factor(c("L-Asparagine", "L-Aspartic acid", "Glycine", "L-Serine", 
               "L-Phenylalanine", "L-Alanine", "L-Lysine", "L-Threonine", 
                "L-Tyrosine", "L-Tryptophan", "L-Proline", "L-Leucine"))

for (metabolite in aminoacyl.metabs.gt){
  mypath<-file.path("aminoacyl", paste("GTboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data<-xc1250b[which(xc1250b$KeggMatch==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", cex.axis=1.5, cex.names=1.5,
          las=2, col=c("grey","grey"),main=paste(metabolite))
  dev.off()
  }
#dev.off()

#box plots for purine metabolism
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots")
setwd(gt_figures)
purine.metabs.gt<-as.factor(c("Urea", "Glycine", "Uric acid", "Inosine", "Adenosine", "Guanosine"))
#pdf("GT_boxplots_purinemetabolism_120315.pdf")
# boxplot for control and 1250 dose for GT metabolites 
# in purine metabolism biosynthesis pathway
for (metabolite in purine.metabs.gt){
  mypath<-file.path("purine", paste("GTboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data<-xc1250b[which(xc1250b$KeggMatch==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", cex.axis=1.5, cex.names=1.5,
          las=2, col=c("grey","grey"),main=paste(metabolite))
  dev.off()
  }
#dev.off()

# box plots for glycine, serine, and threonine metabolism
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots")
setwd(gt_figures)
#pdf("GT_boxplots_glycinemetabolism_120315.pdf")
glycine.metabs.GT<-as.factor(c("L-Tryptophan", "L-Serine", "Glycine", "L-Threonine", "L-Aspartic acid"))
# boxplot for control and 1250 dose for GT metabolites 
# in glycine metabolism biosynthesis pathway
for (metabolite in glycine.metabs.GT){
  mypath<-file.path("glycine", paste("GTboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data<-xc1250b[which(xc1250b$KeggMatch==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", cex.axis=1.5, cex.names=1.5,las=2, 
          col=c("grey","grey"), main=paste(metabolite))
  dev.off()
}
#dev.off()

#box plots for purine/prymidine/arginine/urea metabolism from paper
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots/")
setwd(gt_figures)
purine2.metabs.gt<-as.factor(c("Uric acid", "Inosine", "Adenosine", "Glycine", 
                              "Guanosine", "Creatinine", "glutamate", "D-Ribose",
                              "L-Alanine", "L-Asparagine", "L-Lysine", "L-Proline", 
                              "Ornithine", "Uridine", "Urea"))
# boxplot for control and 1250 dose for GT metabolites 
# in purine metabolism biosynthesis pathway
for (metabolite in purine2.metabs.gt){
  mypath<-file.path("purineagain", paste("GTboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data<-xc1250b[which(xc1250b$KeggMatch==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", cex.axis=1.5, cex.names=1.5,
          las=2, col=c("grey","grey"), main=paste(metabolite))
  dev.off()
}


# create pdf of boxplots in purine/arginine/urea pathway to compare with paper
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots")
setwd(gt_figures)
pdf("GT_boxplots_purine2121715.pdf")
purine2.metabs.gt<-as.factor(c("Glycine", "Uric acid", "Inosine", "Adenosine", 
                               "Guanosine", "Creatinine", "glutamate", "D-Ribose",
                               "L-Alanine", "L-Asparagine", "L-Lysine", "L-Proline", 
                               "Ornithine", "Uridine", "Urea"))
for (metabolite in purine2.metabs.gt){
  temp.data<-xc1250b[which(xc1250b$KeggMatch==metabolite),]
  boxplot(temp.data$value~temp.data$"class$V1", main=paste(metabolite))
}
dev.off()
####################### testing

#################################################################
#### gather data AT ####
#AT_bin_ID_sim <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/AT_bin_ID_sim.csv")
#View(AT_bin_ID_sim)
#AT200bin_IDFinal_sim_1201315 <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/AT200bin_IDFinal_sim_1201315.csv")
#names(AT200bin_IDFinal_sim_1201315)
SVMrank4ID_atp1026154MattID <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/SVMrank4ID_atp1026154MattID.csv")
View(SVMrank4ID_atp1026154MattID)
# import pre-processed bin data 
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/ProcessedTable.csv")
data.at<-cbind(class$V1, ProcessedTable)
#subset data by the 200 top ranked bins with SVM-RFE
data.at.sub <- data.at[,names(data.at) %in% SVMrank4ID_atp1026154MattID$peaks3 ]
dim(data.at.sub)
# add class and name column back
head(names(data.at.sub))
data.at.sub2<-cbind(data.at$X,data.at$`class$V1`,  data.at.sub)
View(data.at.sub2)
# transform data from wide form to long form using lib(reshape2)
names(data.at.sub2)
data.at.r<-melt(data.at.sub2, id.vars=c("data.at$`class$V1`", "data.at$X"))
View(data.at.r)
#rename retention time mass fragment bins with metabolite ID with merge
data.at.r2<-merge(data.at.r,SVMrank4ID_atp1026154MattID, by.x="variable", by.y="peaks3" )
View(data.at.r2)
dim(data.at.r2)
#test<-data.at.r2[which(data.at.r2$MattID1=="l-proline"),]

#drop bins with similiarity value <700
y3<-data.at.r2[data.at.r2$SimMatt1>699,]
dim(y3)
View(y3)
#rename column
names(y3)[2]<-"class"
names(y3)

#########################################################################
#### boxplots AT ####
# box plots for all doses TOF ID
unique.metabolites<-unique(y3$MattID1)
for (metabolite in unique.metabolites){
  temp.data<-y3[which(y3$MattID1==metabolite),]
  boxplot(temp.data$value~temp.data$class, main=paste(metabolite))
}

#subset dataset to control and highest dose
y4<-subset(y3,y3$class =="1" | y3$class=="5" )
View(y4)
dim(y4)
# boxplot for control and 1250 dose TOF ID
unique.metabolites<-unique(y4$MattID1)
for (metabolite in unique.metabolites){
  temp.data<-y4[which(y4$MattID1==metabolite),]
  boxplot(temp.data$value~temp.data$class, main=paste(metabolite))
}
#test<-y4[which(y4$MattID1=="l-proline"),]

#rename metabolites based on KEGG names
#import csv for name translation
AT_Name_Kegg_match120215 <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/AT_Name_Kegg_match120215.csv")
View(AT_Name_Kegg_match120215)
y5<-merge(y4, AT_Name_Kegg_match120215, by.x = "MattID1", by.y = "Query" )
#View(y5)
#print pdf of box plots AT 
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
setwd(gt_figures)
pdf("AT_boxplots.pdf")
# boxplot for control and 1250 dose Kegg ID
unique.metabolites<-unique(y5$Match)
for (metabolite in unique.metabolites){
  temp.data2<-y5[which(y5$Match==metabolite),]
  boxplot(temp.data2$value~temp.data2$class, main=paste(metabolite), is.null=TRUE)
}
dev.off()

####

######################################################################################
#### plot individual box plots as jpg for 3 top pathways from overlap between AT and GT
# box plots for aminoacyl tRNA biosynthesis AT
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots")
setwd(gt_figures)
#pdf("AT_boxplots_aminoacyl_120315.pdf")
# boxplot for control and 1250 dose for AT metabolites 
# in Aminoacyl tRNA biosynthesis pathway
aminoacyl.metabs.AT<-factor(c("Glycine", "L-Serine", "L-Threonine", "L-Proline"))
for (metabolite in aminoacyl.metabs.AT){
  mypath<-file.path("aminoacyl", paste("ATboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data2<-y5[which(y5$Match==metabolite),]
  boxplot(temp.data2$value~temp.data2$class, main=paste(metabolite),cex.axis=1.5, cex.names=1.5,
          las=2, is.null=TRUE)
  dev.off()
  }
#dev.off()

# box plots for purine metabolism AT
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots")
setwd(gt_figures)
#pdf("AT_boxplots_purinemetabolism_120315.pdf")
# boxplot for control and 1250 dose for AT metabolites 
# in Aminoacyl tRNA biosynthesis pathway
purine.metabs.AT<-as.factor(c("Guanosine", "Inosine", "adenosine", "Glycine"))
for (metabolite in purine.metabs.AT){
  mypath<-file.path("purine", paste("ATboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data2<-y5[which(y5$Match==metabolite),]
  boxplot(temp.data2$value~temp.data2$class, main=paste(metabolite), cex.axis=1.5, cex.names=1.5,
          las=2, is.null=TRUE)
  dev.off()
  }
#dev.off()

# box plots for glycine metabolism AT
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots/")
setwd(gt_figures)
#pdf("AT_boxplots_glycinemetabolism_120315.pdf")
# boxplot for control and 1250 dose for AT metabolites 
# in Aminoacyl tRNA biosynthesis pathway
glycine.metabs.AT<-as.factor(c("L-Threonine", "Glycine", "L-Serine"))
for (metabolite in glycine.metabs.AT){
  mypath<-file.path("glycine", paste("ATboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data2<-y5[which(y5$Match==metabolite),]
  boxplot(temp.data2$value~temp.data2$class, main=paste(metabolite), 
          names=c("Control", expression("1250" * ~mu~"g/L")), las=TRUE, cex.axis=1.5, cex.names=1.5,
          is.null=TRUE)
  dev.off()
  }
#dev.off()

##### box plots for purine/arginine/urea metabolism from paper ####
purine2.metabs.at<-as.factor(c("Inosine",  "Glycine", "adenosine",
                               "Guanosine", "Creatinine", "glutamate", "D-Ribose",
                                "L-Proline"))
# boxplot for control and 1250 dose for AT metabolites 
# in purine metabolism biosynthesis pathway
#setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures/boxplots/")
setwd(gt_figures)
for (metabolite in purine2.metabs.at){
  mypath<-file.path("purineagain2", paste("ATboxplot_",metabolite, ".jpg", sep=""))
  jpeg(file=mypath)
  temp.data2<-y5[which(y5$Match==metabolite),]
  boxplot(temp.data2$value~temp.data2$class, cex.axis=1.5, cex.names=1.5,
          las=2, col=c("grey","grey"), main=paste(metabolite))
  dev.off()
}

temp.data2<-y5[which(y5$Match=="adenosine"),] 
temp.data2<-y5[which(y5$Match=="D-Ribose"),] 




# use loop to select metabolite, make box plot of that metabolite
