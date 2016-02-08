library(e1071)
library(muma)
# GT using bins from 1250 biomarker comparison to get classificaiton accuracies 
# from lower doses
# 16November2015 
#### import data normalized with MUMA ####
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz2/WorkDir_allmz9000_noutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz2/WorkDir_allmz9000_noutlier/Preprocessing_Data_a/ProcessedTable.csv")
data2<-cbind(class$V1, ProcessedTable)
View(data2)
#import biomarker bin IDs for subset with GT
biomarkers_gtp041615_MattfinalID2 <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/biomarkers_gtp041615_MattfinalID2.csv")
GT_R4_biomarkersubset <- biomarkers_gtp041615_MattfinalID2
#View(GT_R4_biomarkersubset)
names(GT_R4_biomarkersubset)
#create dataframe from imported data
peaks4mergeGT<-data.frame(GT_R4_biomarkersubset$peaks3, GT_R4_biomarkersubset$ID1)
head(peaks4mergeGT)
summary(peaks4mergeGT)
names(peaks4mergeGT)
#subset whole GT rt mz normalized and scaled dataset by biomarkers identified
x.subGT <- data2[,names(data2) %in% peaks4mergeGT$GT_R4_biomarkersubset.peaks3 ]
View(x.subGT)

#### set up data for c vs. 1250 classification test with biomarkers only ####
x.subGT1250<-x.subGT[-6:-10,]
View(x.subGT1250)
x.subGT1250<-x.subGT1250[-12:-20,]
View(x.subGT1250)
Y<-c(1,1,1,1,1,5,5,5,5,5,5)
#### try tuning 1250 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #89 [5,2]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]
[1,] 38.03636 40.00000 40.00000 35.09091 35.23636 34.14545
[2,] 38.10909 39.45455 38.58182 36.10909 36.03636 35.52727
[3,] 38.50909 39.92727 39.56364 34.94545 34.83636 35.05455
[4,] 35.89091 39.30909 36.69091 33.60000 34.43636 34.58182
[5,] 38.50909 89.01818 36.47273 36.47273 35.05455 34.76364
[6,] 69.30909 88.00000 36.36364 34.87273 34.14545 34.29091
#### try finer tuning 1250 dose model 1st #### 
cost.vector<-seq(1e4,1e6, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #90 [8,10]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
[1,] 86.14545 89.27273 88.98182 88.72727 88.25455 87.70909 88.65455 88.58182
[2,] 88.14545 87.81818 87.85455 87.89091 89.20000 88.00000 88.94545 88.72727
[3,] 88.36364 88.18182 88.25455 88.65455 87.85455 87.30909 87.60000 88.80000
[4,] 88.18182 89.05455 88.58182 88.25455 87.60000 88.72727 89.23636 89.16364
[5,] 89.20000 88.25455 88.72727 88.00000 88.10909 86.43636 88.03636 87.41818
[6,] 87.56364 88.72727 88.36364 87.74545 89.12727 88.21818 88.36364 88.47273
[7,] 88.43636 87.60000 88.90909 88.29091 88.94545 88.69091 88.50909 88.03636
[8,] 88.00000 88.21818 87.85455 89.09091 89.30909 88.36364 89.60000 88.58182
[9,] 88.54545 88.47273 90.03636 88.29091 88.94545 89.45455 87.81818 89.20000
[10,] 88.87273 87.12727 88.10909 88.50909 89.09091 89.23636 88.69091 89.16364
[,9]    [,10]
[1,] 89.09091 87.81818
[2,] 88.14545 87.30909
[3,] 88.18182 89.41818
[4,] 88.36364 88.76364
[5,] 87.81818 88.72727
[6,] 88.76364 88.25455
[7,] 89.05455 88.54545
[8,] 88.61818 90.18182
[9,] 88.14545 87.56364
[10,] 88.72727 88.87273
#cost.value = 780000
#gamma.value = 1.0e-04
#### SVM with dose 1250 repeat 3-fold cross validation after second tuning #### 
accuracy.vectorGT1250<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subGT1250, as.factor(Y), cost = 780000, gamma=1.0e-04, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT1250[i]<-svmModel5$tot.accuracy
}
acc.out.gt.1250<-mean(accuracy.vectorGT1250) #89.34545 % accuracy

#### set up data for c vs. 250 classification test with biomarkers only ####
x.subGT250<-x.subGT[-6:-16,]
View(x.subGT250)
x.subGT250<-x.subGT250[-10:-14,]
View(x.subGT250)
Y<-c(1,1,1,1,1,4,4,4,4) #4 is 250 dose
#### try tuning 250 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #93.1 [6,2]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]
[1,] 53.73333 59.24444 42.44444 43.68889 35.77778 35.60000
[2,] 45.28889 59.46667 46.71111 41.60000 35.91111 36.17778
[3,] 49.60000 57.82222 44.44444 42.48889 36.04444 34.75556
[4,] 48.71111 54.97778 43.15556 47.64444 36.17778 36.35556
[5,] 48.53333 92.22222 45.28889 45.55556 35.28889 35.42222
[6,] 76.04444 93.06667 43.86667 47.15556 36.00000 35.11111
#### try finer tuning 250 dose model 1st #### 
cost.vector<-seq(1e9,1e11,length.out=10)
gamma.vector<-seq(1e-6,1e-4,length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #93.91 [5,9]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
[1,] 93.37778 92.97778 92.00000 91.86667 93.06667 92.26667 92.84444 93.60000
[2,] 93.11111 92.08889 92.66667 92.97778 93.64444 92.62222 92.31111 93.24444
[3,] 92.48889 92.48889 93.55556 92.08889 92.97778 92.84444 92.53333 92.84444
[4,] 92.53333 93.37778 93.46667 92.57778 92.44444 94.08889 92.31111 92.35556
[5,] 92.71111 93.20000 92.35556 93.06667 92.62222 93.06667 92.44444 93.91111
[6,] 93.06667 92.44444 93.37778 92.57778 93.73333 92.84444 93.02222 93.11111
[7,] 94.08889 92.31111 93.20000 92.71111 93.11111 93.11111 92.00000 93.06667
[8,] 92.40000 93.68889 93.15556 93.24444 92.80000 92.48889 93.06667 92.04444
[9,] 93.28889 92.84444 92.40000 93.42222 92.71111 92.26667 93.68889 92.75556
[10,] 93.51111 92.53333 92.84444 93.24444 92.57778 93.37778 93.42222 91.60000
[,9]    [,10]
[1,] 93.82222 93.11111
[2,] 93.02222 92.31111
[3,] 93.95556 92.53333
[4,] 93.86667 92.66667
[5,] 93.91111 92.44444
[6,] 91.91111 93.28889
[7,] 93.73333 91.11111
[8,] 93.33333 93.28889
[9,] 93.11111 93.06667
[10,] 93.37778 92.97778
#cost.value = 4.5e+10
#gamma.value = 8.9e-05
#### SVM with dose 250 repeat 3-fold cross validation after second tuning #### 
accuracy.vectorGT250<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subGT250, as.factor(Y), cost = 4.5e+10, gamma=8.9e-05, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT250[i]<-svmModel5$tot.accuracy
}
acc.out.gt.250<-mean(accuracy.vectorGT250) #92.44444 % accuracy

#### set up data for c vs. 50 classification test with biomarkers only ####
View(data2)
x.subGT50<-x.subGT[-6:-20,]
View(x.subGT50)
Y<-c(1,1,1,1,1,3,3,3,3,3) #4 is 250 dose
#### try tuning 50 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #62.32 [5,3]
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 31.32 33.24 34.88 31.72 31.12 30.76
[2,] 31.92 34.80 34.36 33.48 30.80 30.76
[3,] 31.80 35.72 34.76 34.36 31.08 31.56
[4,] 33.52 34.04 34.16 34.00 31.36 30.80
[5,] 32.92 62.32 33.24 34.16 30.88 31.08
[6,] 59.92 60.92 33.04 34.32 31.16 31.64

#### try finer tuning 50 dose model 1st #### 
cost.vector<-seq(1e4,1e6, length.out=10)
gamma.vector<-seq(1e-2,1,length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #51.16 [4,1]
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
[1,] 50.96 34.28 34.24 33.40 34.20 33.64 34.84 32.40 32.84 34.24
[2,] 47.52 34.44 34.76 33.44 32.76 33.00 33.40 32.84 33.36 33.40
[3,] 48.92 34.48 35.44 35.12 34.44 34.68 34.80 32.80 34.04 33.36
[4,] 51.16 33.88 34.88 34.16 35.52 33.32 34.68 33.48 34.52 33.20
[5,] 49.52 33.40 34.28 34.04 33.72 32.24 33.56 34.32 33.96 33.04
[6,] 49.12 34.60 34.32 34.64 33.80 33.80 34.36 31.72 32.88 33.20
[7,] 50.88 34.88 34.92 34.04 33.76 33.00 35.56 32.76 32.28 33.44
[8,] 50.64 33.40 34.16 33.80 34.72 33.88 32.92 33.68 31.92 32.96
[9,] 49.68 34.20 34.44 33.80 33.72 32.08 33.68 33.00 31.20 32.72
[10,] 49.16 34.40 33.60 34.16 33.64 33.52 34.28 33.84 32.64 33.64
#cost.value = 340000
#gamma.value = 0.01
#### SVM with dose 50 repeat 3-fold cross validation after second tuning #### 
accuracy.vectorGT50<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subGT50, as.factor(Y), cost = 340000, gamma=0.01, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT50[i]<-svmModel5$tot.accuracy
}
acc.out.gt.50<-mean(accuracy.vectorGT50) # 49.28% accuracy

#### set up data for c vs. 10 classification test with biomarkers only ####
View(data2)
x.subGT10<-x.subGT[-11:-25,]
View(x.subGT10)
Y<-c(1,1,1,1,1,2,2,2,2,2) #4 is 250 dose
#### try tuning 10 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #70.68 [6,2]
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 34.56 33.76 32.36 32.80 30.76 30.68
[2,] 33.24 33.88 32.96 33.24 30.24 31.24
[3,] 34.80 32.04 33.64 32.56 30.32 30.92
[4,] 33.72 33.28 33.80 33.44 29.96 30.76
[5,] 32.04 67.04 33.80 32.48 31.80 30.28
[6,] 69.32 70.68 33.04 31.04 30.36 30.04
#### try tuning 10 dose model 1st #### 
cost.vector<-seq(1e11,1e9,length.out=10)
gamma.vector<-seq(1e-6,1e-4,length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subGT10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #70.56 [4,8]
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
[1,] 70.00 69.16 68.04 69.60 68.32 68.48 67.92 67.96 67.00 69.56
[2,] 69.16 69.16 67.92 68.64 70.52 68.40 69.92 68.92 69.32 68.96
[3,] 69.12 67.88 68.48 68.92 68.60 68.92 70.08 67.84 69.32 69.96
[4,] 68.24 68.12 69.16 67.68 68.32 66.80 68.32 70.56 68.60 66.92
[5,] 68.20 68.80 68.80 69.24 68.40 68.32 66.92 68.64 67.96 68.32
[6,] 69.16 70.16 67.00 67.84 68.92 68.96 69.24 68.76 69.08 69.40
[7,] 67.60 70.16 68.84 67.68 68.12 67.76 69.76 68.04 68.24 68.40
[8,] 68.56 67.88 67.56 67.80 66.96 69.64 68.04 69.48 69.96 68.52
[9,] 69.08 69.20 66.48 69.48 68.44 65.28 68.12 68.20 67.24 69.04
[10,] 70.40 68.88 70.36 69.84 69.36 68.32 68.96 69.92 69.24 69.08
#cost.value = 6.7e+10
# gamma.value = 7.8e-05
#### SVM with dose 10 repeat 3-fold cross validation after second tuning #### 
accuracy.vectorGT10<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subGT10, as.factor(Y), cost = 6.7e+10, gamma=7.8e-05, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT10[i]<-svmModel5$tot.accuracy
}
acc.out.gt.10<-mean(accuracy.vectorGT10) # 68.96% accuracy

#### plot of classification accuracy results for biomarker subset with final ID ####
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
pdf("GT_biomarker_accuracy111615.pdf",width=4,height=3)
acc.vectorGT<-c(acc.out.gt.10, acc.out.gt.50, acc.out.gt.250, acc.out.gt.1250)
x.vector<-c(1,2,3,4)
barplot(acc.vectorGT, ylim=c(0,100), yaxt='n', pch=16, xaxt="n", xlab="Dose (ug/L)",ylab="% accuracy")
axis(side = 1, at = seq(1, 4, by = 1), labels = c(10,50,250,1250),  tcl = -0.2)
axis(side = 2, at = seq(0, 100, by = 10), labels = c(0,10,20,30,40,50,60,70,80,90,100),  tcl = -0.2)
dev.off()

########################################################################
#### Look at classification accuracy with SVM-RFE ranked 200 bins across all doses #### 
#### set up data #### 
# import 200 bins data
SVMrank4ID_gtp041615_MattfinalID <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/SVMrank4ID_gtp041615_MattfinalID.csv")
View(SVMrank4ID_gtp041615_MattfinalID)
names(SVMrank4ID_gtp041615_MattfinalID)
#create dataframe from imported data
peaks4mergeGT2<-data.frame(SVMrank4ID_gtp041615_MattfinalID$peaks3, SVMrank4ID_gtp041615_MattfinalID$ID1)
head(peaks4mergeGT2)
summary(peaks4mergeGT2)
names(peaks4mergeGT2)
#subset whole AT rt mz normalized and scaled dataset by biomarkers identified
x.subGT2 <- data2[,names(data2) %in% peaks4mergeGT2$SVMrank4ID_gtp041615_MattfinalID.peaks3 ]
View(x.subGT2)
#### set up data for c vs. 1250 classification test with biomarkers only ####
x.sub2GT1250<-x.subGT2[-6:-10,]
View(x.sub2GT1250)
x.sub2GT1250<-x.sub2GT1250[-12:-20,]
View(x.sub2GT1250)
Y<-c(1,1,1,1,1,5,5,5,5,5,5)
#### try tuning 1250 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
#### output with 2nd tuning for c vs. 1250 #### 
cost.vector<-seq(1e9,1e11, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9

#### SVM with dose 1250 repeat 3-fold cross validation after second tuning #### 
accuracy.vectorGT1250<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.sub2GT1250, as.factor(Y), cost = 2.3e10, gamma=2.3e-5, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT1250[i]<-svmModel5$tot.accuracy
}
acc.out.gt.1250<-mean(accuracy.vectorGT1250) # 88.4% accuracy
#### set up data for c vs. 250 classification test with biomarkers only ####
x.sub2GT250<-x.subGT[-6:-16,]
View(x.sub2GT250)
x.sub2GT250<-x.subGT250[-10:-14,]
View(x.sub2GT250)
Y<-c(1,1,1,1,1,4,4,4,4) #4 is 250 dose

#### try tuning 250 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]
[1,] 43.77778 44.35556 41.15556 36.40000 34.80000 35.33333
[2,] 43.28889 46.00000 41.55556 37.20000 36.08889 34.97778
[3,] 42.17778 44.57778 40.80000 38.40000 35.55556 36.84444
[4,] 40.44444 41.91111 39.91111 38.26667 36.26667 35.51111
[5,] 41.33333 69.77778 40.31111 37.24444 33.46667 35.06667
[6,] 56.71111 69.24444 39.55556 36.80000 34.80000 35.33333

#### try tuning 250 dose model 2nd #### 
cost.vector<-seq(1e4,1e6, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
#### SVM with dose 250 repeat 3-fold cross validation after second tuning #### 
accuracy.vectorGT250<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.sub2GT250, as.factor(Y), cost = 340000, gamma=1e-4, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT250[i]<-svmModel5$tot.accuracy
}
acc.out.gt.250<-mean(accuracy.vectorGT250) # 70% accuracy

#### set up data for c vs. 50 classification test with biomarkers only ####
View(data2)
x.sub2GT50<-x.subGT2[-6:-20,]
View(x.sub2GT50)
Y<-c(1,1,1,1,1,3,3,3,3,3) #3 is 50 dose
#### try tuning 50 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 38.64 37.04 35.92 32.04 32.12 30.84
[2,] 34.48 37.00 34.00 31.92 31.00 30.00
[3,] 34.24 37.04 33.40 30.96 30.68 31.40
[4,] 34.28 36.72 33.72 30.48 32.24 30.88
[5,] 34.92 78.00 31.84 30.40 31.12 29.72
[6,] 67.04 78.24 33.92 31.00 31.60 31.60

#### try tuning 50 dose model 2nd #### 
cost.vector<-seq(1e9,1e11, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9

#### accuracy of 50 vs. control with 200 bins from 1250 ####
accuracy.vectorGT50<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.sub2GT50, as.factor(Y), cost = 5.6e10, gamma=2.3e-5, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT50[i]<-svmModel5$tot.accuracy
}
acc.out.gt.50<-mean(accuracy.vectorGT50) # %77.68

#### set up data for c vs. 10 classification test with top 200 bins from 1250 ####
View(data2)
x.sub2GT10<-x.subGT2[-11:-25,]
View(x.sub2GT10)
Y<-c(1,1,1,1,1,2,2,2,2,2) #2 is 10 dose
#### try tuning 10 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9

#### try tuning 10 dose model 2nd #### 
cost.vector<-seq(1e9,1e11, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub2GT10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9

#### accuracy of 10 vs. control with 200 bins from 1250 ####
accuracy.vectorGT10<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.sub2GT10, as.factor(Y), cost = 1e9, gamma=3.4e-05, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorGT10[i]<-svmModel5$tot.accuracy
}
acc.out.gt.10<-mean(accuracy.vectorGT50) # %77.68 

#### plot of classification accuracy using top 200 bin from 1250 vs. control #### 
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
pdf("GT_200bin1250_accuracy.pdf",width=4,height=3)
acc.vectorGT<-c(acc.out.gt.10, acc.out.gt.50, acc.out.gt.250, acc.out.gt.1250)
x.vector<-c(1,2,3,4)
barplot(acc.vectorGT, ylim=c(0,100), yaxt='n', pch=16, xaxt="n", xlab="Dose (ug/L)",ylab="% accuracy")
axis(side = 1, at = seq(1, 4, by = 1), labels = c(10,50,250,1250),  tcl = -0.2)
axis(side = 2, at = seq(0, 100, by = 10), labels = c(0,10,20,30,40,50,60,70,80,90,100),  tcl = -0.2)
dev.off()

