library(e1071)
library(muma)
# AT using bins from 1250 biomarker comparison to get classificaiton accuracies 
# from lower doses

#### data import for SVM ####
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/ProcessedTable.csv")
data.at<-cbind(class$V1, ProcessedTable)
View(data.at)
dim(data.at)
#data.at$Class<-c(1,1,1,1,5,5,5,5)
# data.at<-data.at[,-1:-2]
# View(data.at)
# dim(data.at)
# Y<-c(1,1,1,1,5,5,5,5)
# Y<-as.factor(Y)

#import biomarker bin IDs for subset with GT
SVM_FinalAT_111615biomarkersubset <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/SVM_FinalAT_111615biomarkersubset.csv")
AT_R4_biomarkersubset<-(SVM_FinalAT_111615biomarkersubset)
View(AT_R4_biomarkersubset)
names(AT_R4_biomarkersubset)
#create dataframe from imported data
peaks4merge<-data.frame(AT_R4_biomarkersubset$peaks3, AT_R4_biomarkersubset$ID1)
head(peaks4merge)
summary(peaks4merge)
names(peaks4merge)
#subset whole AT rt mz normalized and scaled dataset by biomarkers identified
x.sub5 <- data.at[,names(data.at) %in% peaks4merge$AT_R4_biomarkersubset.peaks3 ]
View(x.sub5)
#### set up dataframe for control vs. 250 comparison with biomarker subset ####
x.sub6<-x.sub5[-5:-13,] #delete other doses
View(x.sub6)
x.sub7<-x.sub6[-10:-13,] #delete other doses
View(x.sub7)
Y<-c(1,1,1,1,4,4,4,4,4) #classes of columns 1=control 4=250 dose
#### run SVM #### 
svmModel = svm(x.sub7, as.factor(Y), cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
#svmModel = svm(data.at[, featureRankedList[1:10]], Y, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
svmModel 
summary(svmModel)
print(svmModel)
#### try tuning 250 model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub7, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 # 80% [5,3]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]
[1,] 50.31111 47.86667 47.02222 36.22222 36.13333 35.28889
[2,] 38.71111 47.33333 44.17778 36.88889 35.15556 35.73333
[3,] 37.51111 50.08889 44.13333 35.86667 34.71111 35.46667
[4,] 38.53333 50.04444 43.42222 36.44444 36.04444 36.22222
[5,] 36.71111 80.17778 44.08889 34.84444 37.24444 35.06667
[6,] 65.51111 79.02222 43.60000 34.48889 35.15556 34.97778


#### try finer tuning model ####
cost.vector<-seq(1e4,1e6, length.out=10)
gamma.vector<-seq(1e-2,1,length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub7, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #61% [3,1]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
[1,] 61.42222 43.46667 44.71111 43.15556 45.15556 42.57778 46.13333 44.00000
[2,] 60.22222 44.62222 43.55556 43.28889 40.88889 43.60000 41.42222 43.33333
[3,] 61.95556 44.88889 43.06667 44.35556 44.80000 42.88889 44.53333 45.46667
[4,] 61.73333 44.35556 45.20000 45.11111 43.15556 45.42222 44.93333 44.08889
[5,] 60.31111 45.55556 44.97778 43.86667 43.82222 41.86667 46.13333 43.86667
[6,] 60.31111 44.40000 43.60000 43.02222 44.62222 45.28889 43.95556 45.77778
[7,] 59.06667 44.08889 44.31111 42.88889 43.33333 42.00000 45.20000 46.04444
[8,] 59.60000 44.44444 46.44444 45.02222 42.13333 42.04444 42.66667 44.35556
[9,] 59.73333 43.77778 45.37778 46.08889 44.13333 42.88889 43.64444 43.33333
[10,] 61.60000 44.88889 44.31111 43.86667 43.24444 43.60000 44.97778 46.53333
[,9]    [,10]
[1,] 45.51111 45.02222
[2,] 42.88889 48.53333
[3,] 44.62222 46.93333
[4,] 44.31111 47.77778
[5,] 47.15556 44.48889
[6,] 46.62222 47.15556
[7,] 43.82222 45.51111
[8,] 44.04444 46.22222
[9,] 47.02222 44.48889
[10,] 45.82222 45.68889
#### SVM with 250 repeat 3-fold cross validation after second tuning #### 
accuracy.vector<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel3 = svm(x.sub7, as.factor(Y), cost = 230000, gamma=0.01, kernel="radial", cross=3, scale=F )
  accuracy.vector[i]<-svmModel3$tot.accuracy
}
acc.out.at.250<-mean(accuracy.vector) #56.88889% accuracy

#### format data to do SVM for biomarker subset on 50 dose ####
View(x.sub5)
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/class.csv")
View(class)
x.sub50<-x.sub5[-5:-18,] #delete other doses
View(x.sub50)
Y<-c(1,1,1,1,3,3,3,3) #classes of columns 1=control 3=50 dose
#### run test SVM 50#### 
svmModel50 = svm(x.sub50, as.factor(Y), cost = 8.9e9, gamma=8.9e-6, kernel="radial", cross=3 )
svmModel50 
summary(svmModel50)
print(svmModel50)
#### try tuning 50 model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 # 54.45% [6,3]
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 27.90 27.95 30.15 27.60 29.05 28.55
[2,] 27.75 28.80 30.05 29.05 28.45 28.40
[3,] 28.50 26.95 29.65 28.40 27.95 29.20
[4,] 28.20 27.45 28.25 29.25 28.45 29.75
[5,] 29.45 51.65 29.10 28.10 28.45 28.60
[6,] 42.25 54.45 28.20 28.45 28.65 28.75

#### try tuning 50 model again finer #### 
cost.vector<-seq(1e9,1e11, length.out=10)
gamma.vector<-seq(1e-2,1, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #43.15 [4,1]
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
[1,] 40.15 28.55 30.35 30.00 30.85 30.80 31.45 29.25 29.60 30.10
[2,] 42.65 30.60 30.30 30.90 29.95 30.70 30.35 29.30 28.85 29.70
[3,] 41.85 27.95 29.70 30.95 29.45 29.55 32.00 29.35 29.30 29.35
[4,] 43.15 30.80 30.75 30.20 30.15 28.10 31.65 30.15 30.35 30.70
[5,] 42.90 29.75 31.10 30.15 30.85 29.45 29.50 30.20 30.95 29.85
[6,] 41.30 30.00 29.70 30.75 28.20 30.65 29.30 30.55 31.00 30.40
[7,] 42.75 28.80 30.30 28.35 30.10 29.15 30.15 31.40 29.55 30.25
[8,] 41.20 30.00 30.65 31.65 30.50 31.35 30.10 30.15 30.90 30.25
[9,] 43.00 29.90 31.45 29.50 30.10 30.10 29.75 30.45 28.30 30.05
[10,] 40.05 29.50 31.50 31.00 29.95 29.45 30.60 29.60 29.90 31.25

#### SVM with 50 repeat 3-fold cross validation after second tuning #### 
accuracy.vector<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel3 = svm(x.sub50, as.factor(Y), cost = 3.4e+10, gamma=0.01, kernel="radial", cross=3 )
  accuracy.vector[i]<-svmModel3$tot.accuracy
}
acc.out.at.50<-mean(accuracy.vector) #35.45 % accuracy

#### format data to do SVM for biomarker subset on 10 dose ####
View(x.sub5)
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/class.csv")
View(class)
x.sub10<-x.sub5[-10:-23,] #delete other doses
dim(x.sub5)
View(x.sub10)
Y<-c(1,1,1,1,2,2,2,2,2) #classes of columns 1=control 2=10 dose

#### try tuning 10 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #[6,2] 47%
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]
[1,] 36.97778 34.04444 37.91111 33.64444 35.15556 34.00000
[2,] 32.93333 34.66667 39.64444 36.08889 36.35556 36.40000
[3,] 33.33333 34.88889 41.02222 34.80000 35.20000 33.95556
[4,] 31.91111 36.48889 41.02222 34.17778 36.62222 37.02222
[5,] 33.02222 46.80000 40.53333 36.66667 36.40000 35.77778
[6,] 45.24444 47.55556 38.04444 34.44444 36.13333 34.53333
#### try finer tuning model dose 10####
cost.vector<-seq(1e9,1e11, length.out=10)
gamma.vector<-seq(1e-6,1e-4,length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #47.60000 [5,1]
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
[1,] 46.26667 47.20000 45.95556 44.75556 46.04444 45.51111 45.42222 45.91111
[2,] 46.84444 45.24444 45.37778 47.55556 45.24444 46.26667 45.73333 45.46667
[3,] 46.93333 44.80000 47.33333 47.28889 44.13333 46.57778 44.75556 45.24444
[4,] 45.06667 44.88889 45.68889 46.04444 44.88889 46.13333 45.33333 45.11111
[5,] 47.60000 45.46667 44.71111 46.80000 46.31111 46.53333 46.66667 44.31111
[6,] 45.91111 45.42222 46.57778 44.35556 45.06667 45.55556 44.71111 45.60000
[7,] 45.86667 47.33333 44.97778 46.31111 44.88889 44.57778 43.91111 44.44444
[8,] 44.53333 45.64444 44.13333 46.00000 45.60000 46.13333 46.22222 45.95556
[9,] 47.15556 44.08889 44.44444 46.62222 46.40000 45.37778 45.64444 44.35556
[10,] 46.93333 46.40000 44.97778 45.24444 46.93333 46.44444 45.28889 44.97778
[,9]    [,10]
[1,] 46.66667 43.28889
[2,] 46.35556 47.11111
[3,] 46.31111 45.37778
[4,] 45.82222 46.44444
[5,] 45.20000 44.97778
[6,] 45.37778 45.46667
[7,] 45.33333 46.26667
[8,] 46.08889 46.66667
[9,] 45.86667 45.91111
[10,] 45.46667 44.22222
#### SVM with 10 dose 250 repeat 3-fold cross validation after second tuning #### 
accuracy.vector<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel4 = svm(x.sub10, as.factor(Y), cost = 4.5e+10, gamma=1.0e-06, kernel="radial", cross=3 )
  accuracy.vector[i]<-svmModel4$tot.accuracy
}
acc.out.at.10<-mean(accuracy.vector) # 44.62222% accuracy

#### format data to do SVM for biomarker subset on 1250 dose ####
View(x.sub5)
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/class.csv")
View(class)
x.sub1250<-x.sub5[-5:-9,] #delete other doses
View(x.sub1250)
x.sub1250<-x.sub1250[-9:-17,]
Y<-c(1,1,1,1,5,5,5,5)
#### try tuning 1250 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #[6,2] 64.9%

[,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 32.55 30.15 31.45 28.60 27.15 27.35
[2,] 29.45 29.10 30.00 28.55 26.80 29.05
[3,] 30.40 31.30 29.25 28.80 28.90 29.35
[4,] 30.55 29.35 29.90 28.60 27.85 28.85
[5,] 31.90 64.95 28.05 28.75 29.00 28.45
[6,] 52.75 64.90 31.10 28.80 28.35 28.75

#### try finer tuning 1250 dose model 1st #### 
cost.vector<-seq(1e9,1e11, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.sub1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #[8,2] 67.15%
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
[1,] 63.95 65.45 63.05 66.20 64.00 60.60 65.85 63.90 64.95 63.10
[2,] 62.60 63.35 62.40 65.75 65.25 64.15 64.95 64.80 64.40 61.70
[3,] 63.30 64.50 63.40 66.05 64.05 64.90 65.25 65.50 65.00 63.95
[4,] 66.15 65.95 64.55 63.00 63.20 64.55 62.65 65.40 61.75 61.95
[5,] 66.40 65.20 64.65 64.30 64.85 64.40 64.05 63.60 63.75 63.35
[6,] 64.50 65.75 65.30 63.80 62.75 66.20 64.85 66.65 64.50 63.60
[7,] 63.30 63.65 67.55 64.60 66.25 65.15 65.25 65.20 63.70 65.05
[8,] 63.00 67.15 63.30 65.55 65.35 63.75 62.20 62.50 65.10 65.70
[9,] 64.15 65.75 64.05 65.95 63.85 65.70 65.00 63.80 64.30 64.05
[10,] 65.40 65.60 64.85 66.55 62.85 64.85 61.95 65.45 63.45 63.80
#### try finer tuning again w/ 1250 dose model 1st #### 
#cost.vector<-seq(8.0e5,9.0e5 length.out=10)
#gamma.vector<-seq(4.6e-5,6.6e-6, length.out=10)
#accuracy.vector<-NULL
#m9<-matrix(NA, length(cost.vector), length(gamma.vector))

#for (k in seq_along(cost.vector)){
#  for (j in seq_along(gamma.vector)){    
#    for (i in 1:250){
#      model.radial.all <- svm(x.sub1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
#      accuracy.vector[i]<-model.radial.all$tot.accuracy
#    }
#    m9[k, j]<-mean(accuracy.vector)  
#  }
#}
#m9
#### SVM with 10 dose 1250 repeat 3-fold cross validation after second tuning #### 
accuracy.vector<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.sub1250, as.factor(Y), cost = 7.8e+10, gamma=1.2e-05, kernel="radial", cross=3 )
  accuracy.vector[i]<-svmModel5$tot.accuracy
}
acc.out.at.1250<-mean(accuracy.vector) # 64.15% accuracy

#### plot of classification accuracy results ####
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
pdf("AT_biomarker_accuracy2.pdf",width=4,height=3)
acc.vectorAT<-c(acc.out.at.10,acc.out.at.50,acc.out.at.250,acc.out.at.1250)
x.vector<-c(1,2,3,4)
barplot(acc.vectorAT, ylim=c(0,100), yaxt='n', pch=16, xaxt="n", xlab="Dose (ug/L)",ylab="% accuracy")
axis(side = 1, at = seq(1, 4, by = 1), labels = c(10,50,250,1250),  tcl = -0.2)
axis(side = 2, at = seq(0, 100, by = 10), labels = c(0,10,20,30,40,50,60,70,80,90,100),  tcl = -0.2)
dev.off()

setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
jpeg('AT_biomarker_accuracy2.jpg', quality = 100, bg = "white", res = 300, width = 7, height = 7, units = "in")
acc.vectorAT<-c(acc.out.at.10,acc.out.at.50,acc.out.at.250,acc.out.at.1250)
x.vector<-c(1,2,3,4)
barplot(acc.vectorAT, ylim=c(0,100), yaxt='n', pch=16, xaxt="n", xlab="Dose (ug/L)",ylab="% accuracy")
axis(side = 1, at = seq(1, 4, by = 1), labels = c(10,50,250,1250),  tcl = -0.2)
axis(side = 2, at = seq(0, 100, by = 10), labels = c(0,10,20,30,40,50,60,70,80,90,100),  tcl = -0.2)
dev.off()

pdf("AT_GT_biomarker_accuracy2.pdf",width=4,height=3)
counts <- rbind(acc.vectorAT, acc.vectorGT)
barplot(counts, ylim=c(0,100), yaxt='n', axis.lty=1, pch=16,  
        xlab="Dose (ug/L)",ylab="% accuracy", beside=TRUE)
axis(side = 1, at = seq(2, 12, by = 3), labels = c(10,50,250,1250),  tcl = -0.2)
axis(side = 2, at = seq(0, 100, by = 10), labels = c(0,10,20,30,40,50,60,70,80,90,100),  tcl = -0.2)
dev.off()

setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
jpeg('AT_GT_biomarker_accuracy2.jpg', quality = 100, bg = "white", res = 300, width = 7, height = 7, units = "in")
counts <- rbind(acc.vectorAT, acc.vectorGT)
barplot(counts, ylim=c(0,100), yaxt='n', axis.lty=1, pch=16,  
        xlab="Dose (ug/L)",ylab="% accuracy", beside=TRUE)
axis(side = 1, at = seq(2, 12, by = 3), labels = c(10,50,250,1250),  tcl = -0.2)
axis(side = 2, at = seq(0, 100, by = 10), labels = c(0,10,20,30,40,50,60,70,80,90,100),  tcl = -0.2)
dev.off()

########################################
#### import data for classification accuracy test with 200 bins for AT ####
#########################################
#### reimport data to make sure have scaled normalized data ####
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/ProcessedTable.csv")
data3<-cbind(class$V1, ProcessedTable)
#import 200 bin IDs for subset with AT
SVMrank200_4class110815 <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/SVMrank200_4class110815.csv")
names(SVMrank200_4class110815)
peaks4subset<-as.data.frame(SVMrank200_4class110815)
#subset whole AT rt mz normalized and scaled dataset by 200 ranked bins identified
x.subAT <- data3[,names(data3) %in% peaks4subset$peaks3 ]
View(x.subAT)
dim(x.subAT) #22 200
dim(data3)
names(data3)
View(data3)
#### set up data for c vs. 1250 classification test with subset only ####
x.subAT1250<-x.subAT[-5:-9,]
View(x.subAT1250)
x.subAT1250<-x.subAT1250[-9:-17,]
View(x.subAT1250)
Y<-c(1,1,1,1,5,5,5,5)
#### try tuning 1250 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #58
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 31.00 32.50 30.25 28.05 28.20 28.20
[2,] 29.10 29.30 30.60 29.60 28.25 29.15
[3,] 30.40 32.25 29.00 28.10 26.75 28.80
[4,] 28.30 30.35 29.60 29.55 28.15 28.85
[5,] 29.55 57.10 29.70 29.40 30.05 27.75
[6,] 58.30 58.80 29.65 27.80 28.65 28.60
#### try tuning 1250 dose model again #### 
cost.vector<-seq(1e9,1e11, length.out=10)
gamma.vector<-seq(1e-11,1e-9, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT1250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #65
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
[1,] 53.15 60.50 51.95 58.60 53.55 60.25 61.60 58.65 58.65 57.25
[2,] 60.05 62.00 53.50 58.55 56.35 57.85 59.35 60.30 60.20 58.05
[3,] 61.75 62.65 53.05 58.30 53.55 59.95 61.50 58.55 60.75 59.00
[4,] 62.45 60.35 52.60 60.10 55.50 59.60 60.60 58.40 58.50 59.25
[5,] 62.00 62.95 51.30 58.10 55.90 58.45 62.15 58.80 59.15 55.45
[6,] 61.80 62.10 51.10 57.85 55.65 57.05 63.65 58.70 59.95 59.65
[7,] 63.40 61.20 51.00 57.40 55.50 57.50 61.20 58.30 57.80 57.00
[8,] 64.85 61.30 52.05 58.00 54.95 60.00 60.50 58.10 58.60 58.15
[9,] 62.85 60.95 53.25 57.75 55.30 59.60 60.45 59.35 59.50 57.85
[10,] 61.15 61.45 51.05 56.25 55.85 57.70 59.45 59.55 56.85 56.60

#### get classification accuracy for 1250 ####
accuracy.vectorAT1250<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subAT1250, as.factor(Y), cost = 7.8e+10, gamma=1.0e-11, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorAT1250[i]<-svmModel5$tot.accuracy
}
acc.out.at.1250<-mean(accuracy.vectorAT1250) #64.05%

#### set up data for c vs. 250 classification test with subset only ####
x.subAT250<-x.subAT[-5:-13,]
View(x.subAT250)
x.subAT250<-x.subAT250[-10:-13,]
View(x.subAT250)
Y<-c(1,1,1,1,4,4,4,4,4)
#### try tuning 250 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9 #73%

[,1]     [,2]     [,3]     [,4]     [,5]     [,6]
[1,] 46.44444 46.22222 49.24444 35.95556 34.97778 35.55556
[2,] 46.40000 46.62222 44.08889 37.46667 35.60000 35.95556
[3,] 45.15556 46.93333 44.04444 36.84444 35.91111 36.08889
[4,] 45.95556 45.51111 42.53333 35.06667 36.57778 35.06667
[5,] 43.82222 72.44444 42.66667 36.84444 33.95556 35.91111
[6,] 60.84444 72.22222 41.28889 36.04444 36.53333 36.48889
#### try tuning 250 dose model 2nd #### 

cost.vector<-seq(1e4,1e6, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT250, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
[1,] 70.66667 71.68889 71.91111 72.35556 70.66667 73.02222 70.40000 71.82222
[2,] 72.66667 69.60000 72.40000 72.44444 70.00000 71.15556 69.37778 70.00000
[3,] 72.44444 70.75556 72.35556 72.84444 72.04444 70.04444 71.15556 71.86667
[4,] 73.15556 70.84444 71.51111 71.46667 70.75556 71.46667 71.73333 69.86667
[5,] 72.22222 72.08889 73.24444 72.13333 72.66667 71.55556 69.91111 71.77778
[6,] 71.64444 69.24444 71.55556 69.91111 71.11111 71.11111 72.48889 70.17778
[7,] 70.57778 70.88889 71.15556 73.68889 70.71111 72.04444 70.53333 70.66667
[8,] 70.22222 70.31111 71.82222 72.04444 72.97778 69.86667 72.88889 71.51111
[9,] 72.97778 71.42222 71.68889 72.93333 69.28889 70.22222 70.40000 72.00000
[10,] 72.04444 72.35556 71.91111 71.51111 70.04444 69.82222 71.95556 71.06667
[,9]    [,10]
[1,] 72.31111 71.15556
[2,] 70.71111 72.53333
[3,] 70.53333 70.80000
[4,] 70.17778 71.68889
[5,] 70.17778 71.64444
[6,] 71.33333 70.71111
[7,] 71.60000 70.88889
[8,] 70.35556 71.37778
[9,] 72.17778 71.73333
[10,] 71.02222 69.11111

#### get classification accuracy for 250 ####
accuracy.vectorAT250<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subAT250, as.factor(Y), cost = 670000, gamma=3.4e-05, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorAT250[i]<-svmModel5$tot.accuracy
}
acc.out.at.250<-mean(accuracy.vectorAT250) #71.6

#### set up data for c vs. 50 classification test with subset only ####
x.subAT50<-x.subAT[-5:-18,]
View(x.subAT50)
Y<-c(1,1,1,1,3,3,3,3)

#### try tuning 50 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]
[1,] 27.75 26.10 31.30 29.35 27.35 28.75
[2,] 27.50 27.55 30.35 28.50 29.10 28.05
[3,] 25.85 28.15 30.20 28.60 29.25 29.20
[4,] 25.50 28.65 30.15 29.15 28.85 28.20
[5,] 25.95 49.60 30.65 28.10 27.85 28.80
[6,] 40.55 47.35 30.05 29.05 28.35 29.35
#### try tuning 50 dose model 2nd #### 
cost.vector<-seq(1e4,1e6, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT50, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
[,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
[1,] 47.75 48.20 46.65 46.20 50.25 44.95 47.90 47.50 46.75 46.00
[2,] 49.85 50.10 49.25 46.75 49.80 47.80 47.50 48.45 48.85 47.15
[3,] 47.25 45.90 48.90 49.15 48.90 48.40 48.20 46.25 47.40 48.65
[4,] 49.20 48.65 49.10 47.10 49.20 50.40 47.30 50.90 48.90 46.00
[5,] 47.35 48.20 48.45 47.70 49.55 49.45 48.85 48.20 46.85 46.00
[6,] 48.55 49.50 49.05 47.65 45.10 46.90 47.80 46.55 43.90 47.50
[7,] 49.05 49.50 47.10 47.00 48.60 46.70 47.05 47.10 45.75 45.80
[8,] 48.90 49.70 46.45 48.30 47.90 49.45 49.65 47.05 46.75 48.75
[9,] 48.60 49.30 47.20 50.35 48.60 49.00 47.10 46.25 47.35 47.15
[10,] 48.65 48.25 44.85 46.75 46.25 45.80 49.45 47.20 45.85 47.10
#### get classification accuracy for 50 ####
accuracy.vectorAT50<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subAT50, as.factor(Y), cost = 890000, gamma=3.4e-05, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorAT50[i]<-svmModel5$tot.accuracy
}
acc.out.at.50<-mean(accuracy.vectorAT50) #48.3

#### set up data for c vs. 10 classification test with subset only ####
x.subAT10<-x.subAT[-10:-22,]
View(x.subAT10)
Y<-c(1,1,1,1,2,2,2,2,2)

#### try tuning 10 dose model 1st #### 
cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]
[1,] 35.55556 32.26667 41.33333 37.33333 36.97778 34.53333
[2,] 37.60000 33.86667 42.80000 35.64444 35.28889 35.42222
[3,] 38.44444 34.31111 41.15556 35.95556 35.46667 35.95556
[4,] 39.33333 34.66667 40.48889 33.33333 36.22222 35.46667
[5,] 37.91111 49.42222 39.02222 36.22222 36.26667 36.26667
[6,] 39.20000 52.08889 40.40000 36.22222 35.02222 34.66667

#### try tuning 10 dose model 2nd #### 
cost.vector<-seq(1e11,1e9, length.out=10)
gamma.vector<-seq(1e-6,1e-4, length.out=10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(x.subAT10, as.factor(Y), kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
[,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
[1,] 49.37778 49.95556 49.37778 50.48889 49.11111 49.77778 49.82222 49.28889
[2,] 49.55556 49.68889 50.48889 50.80000 50.88889 51.24444 50.00000 50.40000
[3,] 50.35556 48.71111 49.37778 51.77778 48.40000 49.20000 49.06667 49.20000
[4,] 49.37778 50.31111 50.53333 49.37778 48.31111 49.24444 49.24444 48.22222
[5,] 49.73333 48.53333 50.17778 48.62222 50.88889 50.57778 49.11111 49.15556
[6,] 48.35556 50.00000 52.04444 50.71111 50.22222 50.75556 51.77778 50.35556
[7,] 48.75556 48.75556 52.13333 49.82222 49.95556 51.60000 50.57778 49.64444
[8,] 51.11111 48.40000 50.53333 50.35556 49.82222 50.31111 51.55556 52.40000
[9,] 51.37778 51.15556 51.28889 49.60000 49.91111 49.77778 50.35556 50.71111
[10,] 50.57778 48.13333 49.60000 48.66667 51.95556 50.26667 49.73333 51.28889
[,9]    [,10]
[1,] 49.51111 49.02222
[2,] 48.66667 50.35556
[3,] 49.28889 47.82222
[4,] 48.13333 49.15556
[5,] 47.02222 49.73333
[6,] 46.97778 47.42222
[7,] 50.84444 48.62222
[8,] 47.64444 49.02222
[9,] 49.42222 48.17778
[10,] 50.75556 48.75556
#### get classification accuracy for 10 ####
accuracy.vectorAT10<-NULL
out.acc<-NULL
for (i in 1:250){
  svmModel5 = svm(x.subAT10, as.factor(Y), cost = 2.3e+10, gamma=7.8e-05, kernel="radial", scale=FALSE, cross=3 )
  accuracy.vectorAT10[i]<-svmModel5$tot.accuracy
}
acc.out.at.10<-mean(accuracy.vectorAT10) #50.48889

#### plot of classification accuracy results ####
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
pdf("AT_200subset_accuracy.pdf",width=4,height=3)
acc.vectorAT<-c(acc.out.at.10, acc.out.at.50, acc.out.at.250, acc.out.at.1250)
x.vector<-c(1,2,3,4)
barplot(acc.vectorAT, ylim=c(0,100), yaxt='n', pch=16, xaxt="n", xlab="Dose (ug/L)",ylab="% accuracy")
axis(side = 1, at = seq(1, 4, by = 1), labels = c(10,50,250,1250),  tcl = -0.2)
axis(side = 2, at = seq(0, 100, by = 10), labels = c(0,10,20,30,40,50,60,70,80,90,100),  tcl = -0.2)
dev.off()
