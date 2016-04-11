library(e1071)
library(muma)
#### data import and format for MUMA #### 
XCMSdiff_at_tof_r4 <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/XCMSdiff_at_tof_r4.csv")
atall<-XCMSdiff_at_tof_r4
n<-atall$name
atall_t<-as.data.frame(t(atall[,-1]))
colnames(atall_t)<-n
View(atall_t)
atall_t$Samples <- factor(colnames(atall)[-1])
[1] ca3_1       cb3__1      cc2__1      cd3__1     
[5] ce1_1       X10_a1_1    X10_b3_1    X10_c3__1  
[9] X10_d1_1    X10_e1_1    X1250_a4__1 X1250_b3__1
[13] X1250_c3__1 X1250_c4_1  X1250_d1_1  X250_a2_1  
[17] X250_a4_1   X250_b3_1   X250_c3__1  X250_d2_1  
[21] X250_e1_1   X50_a1_1    X50_b3__1   X50_d3__1  
[25] X50_e1_1   
atall_t$Class<-c("1","1","1","1","1", "2", "2","2","2","2","5","5","5","5","5",
               "4","4","4","4","4","4","3","3","3","3")
View(atall_t)
head(summary(atall_t))
atall_t<-atall_t[sapply(atall_t, function(atall_t) !any(is.na(atall_t)))]#get rid of columns with NA
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615")
write.csv(atall_t, file="atall_t_out.csv")

#### MUMA pca for outliers ####
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615")
work.dir(dir.name="WorkDir_allmz")
# delete first row in excel of .csv and create Class column
explore.data(file="atall_t_out.csv", scaling="a", scal=TRUE, normalize=TRUE, imputation=FALSE, imput="ImputType")
par( mfrow = c( 1, 2 ) )
Plot.pca(pcx=1, pcy=2, scaling="a", test.outlier=TRUE)
Plot.pca(pcx=1, pcy=3, scaling="a", test.outlier=TRUE)
explore.data(file="atall_t_out.csv", scaling="p", scal=TRUE, normalize=TRUE, imputation=FALSE, imput="ImputType")
par( mfrow = c( 1, 2 ) )
Plot.pca(pcx=1, pcy=2, scaling="p", test.outlier=TRUE)
Plot.pca(pcx=1, pcy=3, scaling="p", test.outlier=TRUE)
# drop outliers
atall_t_out <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/atall_t_out.csv")
atall_t_out$Samples
atall_t_out<-atall_t_out[-5,] #ce1_1
atall_t_out<-atall_t_out[-10,] #250c3

write.csv(atall_t_out, file="atall_t2_out.csv")

#### PCA and PLS-da without ce1_1 #### 
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615")
work.dir(dir.name="WorkDir_allmz_wooutlier")
explore.data(file="atall_t2_out.csv", scaling="a", scal=TRUE, normalize=TRUE, imputation=FALSE, imput="ImputType")
par( mfrow = c( 1, 2 ) )
Plot.pca(pcx=1, pcy=2, scaling="a", test.outlier=TRUE)
Plot.pca(pcx=1, pcy=3, scaling="a", test.outlier=TRUE)
plsda(scaling="a") 
Plot.plsda(pcx=1, pcy=2, scaling="a")
Plot.plsda(pcx=2, pcy=3, scaling="a")
explore.data(file="atall_t2_out.csv", scaling="p", scal=TRUE, normalize=TRUE, imputation=FALSE, imput="ImputType")
par( mfrow = c( 1, 2 ) )
Plot.pca(pcx=1, pcy=2, scaling="p", test.outlier=TRUE)
Plot.pca(pcx=1, pcy=3, scaling="p", test.outlier=TRUE)
plsda(scaling="p") 
Plot.plsda(pcx=1, pcy=2, scaling="p")

#### 111815 to create figure of SVM accuracy by number of bins included ####
#### data import for SVM ####
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/WorkDir_allmz_wooutlier/Preprocessing_Data_a/ProcessedTable.csv")

data.at<-cbind(class$V1, ProcessedTable)
View(data.at)
data.at<-data.at[-5:-9,]
View(data.at)
data.at<-data.at[-9:-17,]
View(data.at)
#data2$Class<-as.factor(data2$class$V1)
data.at$Class<-c(1,1,1,1,5,5,5,5)
data.at<-data.at[,-1:-2]
View(data.at)
dim(data.at)
Y<-c(1,1,1,1,5,5,5,5)
Y<-as.factor(Y)
#### other R implementation of RFE for SVM ####
svmrfeFeatureRanking = function(x,y){
  n = ncol(x)
  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 6.4e9, gamma=1.0e-10, cachesize=500,
                   scale=F, type="C-classification", kernel="radial" )
    #compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV
    #compute ranking criteria
    rankingCriteria = w * w
    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix
    #update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
  }
  return (featureRankedList)
}

#### SVM feature ranking #### 
featureRankedList <-svmrfeFeatureRanking(data.at,Y)
# test
svmModel = svm(data.at[, featureRankedList[1:10]], Y, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
svmModel 
summary(svmModel)
dim(data.at)
#### for accuracy plot AT ####
accuracy.vector.at<-NULL
no.features.at<-seq(1,11700,by=100)
out.acc.at<-NULL
for (j in seq_along(no.features.at)){
  for (i in 1:250){
    svmModel.at = svm(data.at[, featureRankedList[1:no.features.at[j]]], Y, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
    accuracy.vector.at[i]<-svmModel.at$tot.accuracy
  }
  out.acc.at[j]<-mean(accuracy.vector.at)  
}
out.acc.at
plot(out.acc.at~no.features.at, xlab="Number of bins", ylab="% Accuracy", las=1, ylim=c(0,100), pch=16)
peaks.at<-t(Y)
View(peaks.at)

setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
png("AT_accuracyperbin.png",height = 4, width = 4, units = 'in', res=600)
plot(out.acc.at~no.features.at, xlab="Number of bins", ylab="% Accuracy", las=1, ylim=c(0,100), pch=16)
dev.off()

#### for second accuracy plot to get a finer look at 1-250 ####
accuracy.vector.at2<-NULL
no.features.at2<-seq(1,10,by=1)
out.acc.at2<-NULL
for (j in seq_along(no.features.at2)){
  for (i in 1:250){
    svmModel.at = svm(data.at[, featureRankedList[1:no.features.at2[j]]], Y, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
    accuracy.vector.at2[i]<-svmModel.at$tot.accuracy
  }
  out.acc.at2[j]<-mean(accuracy.vector.at2)  
}
out.acc.at2
plot(out.acc.at2~no.features.at2, xlab="Number of bins", ylab="% Accuracy", las=1, ylim=c(0,100), pch=16)


#### ttest to compare p values with rfe ####
XCMSdiff_at_tof_r4 <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/MUMA_allpeaks061615/XCMSdiff_at_tof_r4.csv")
View(XCMSdiff_at_tof_r4)
dim(XCMSdiff_at_tof_r4)

ttest_at<-XCMSdiff_at_tof_r4
dim(ttest_at)
# delete unnecessary columns
names(ttest_at)
ttest_at<-ttest_at[,-1]
names(ttest_at)
dim(ttest_at)
names(ttest_at)
ttest_at2<-ttest_at[,-6:-10]
names(ttest_at2)
# drop 10,50,250 doses
ttest_at2<-ttest_at2[,-11:-20]
names(ttest_at2)
# get peak names
XCMSdiff_at_tof_r4_mzrt <- read.csv("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold/XCMSdiff_at_tof_r4_mzrt.csv")
View(XCMSdiff_at_tof_r4_mzrt)
peaks_at<-XCMSdiff_at_tof_r4_mzrt$name
length(peaks_at)
#
# normalize each column sum column abundance/total
ttest.norm.at<-matrix(NA, nrow(ttest_at2), ncol(ttest_at2))
for (i in 1:ncol(ttest_at2)){
  outsum2<-sum(ttest_at2[,i])
  for (j in 1:nrow(ttest_at2)){
    ttest.norm.at[j,i]<-ttest_at2[j,i] / outsum2
  }
}
View(ttest.norm.at)
ttest.norm.at2<-cbind(peaks_at,as.data.frame(ttest.norm.at))
View(ttest.norm.at2)
dim(ttest.norm.at2)
names(ttest_at2)
names2<-c("rtmz","ca3_1","cb3__1","cc2__1","cd3__1","ce1_1","X1250_a4__1",
        "X1250_b3__1","X1250_c3__1","X1250_c4_1","X1250_d1_1")
colnames(ttest.norm.at2)<-names2
View(ttest.norm.at2)

#transform matrix 
ttest.norm.at3<-t(ttest.norm.at2)
View(ttest.norm.at3)
colnames(ttest.norm.at3)<-ttest.norm.at3[1,]
View(ttest.norm.at3)
ttest.norm.at3<-ttest.norm.at3[-1,]
View(ttest.norm.at3)
Class<-c("c","c","c","c","c","d1250","d1250","d1250","d1250","d1250")
Class<-as.factor(Class)
class(Class)
ttest.norm.at4<-cbind(Class, ttest.norm.at3)
View(ttest.norm.at4)
ttest.norm.at4<-as.data.frame(ttest.norm.at4)
View(ttest.norm.at4)
class(ttest.norm.at4)

#convert to numeric
test<-as.numeric(as.character(ttest.norm.at4[,3]))
class(test)
View(test)
class(ttest.norm.at4[,5])
ttest.norm.at4<-as.numeric((unlist(ttest.norm.at4[,-1]))                 
                        for( i in 2:ncol(ttest.norm.at4)){
                          ttest.norm.at4[,i]<-as.numeric(as.character(ttest.norm.at4[,i]))
                        }
summary(ttest.norm.at4[,5])
                           
# loop through tttest of each row
t.out.test<-t.test(ttest.norm.at4[,6]~ttest.norm.at4$Class )
summary(t.out.test)                        
t.out.test$p.value
#set up empty vector
dim(ttest.norm.at4)
p.out.at<-vector()
#loop through each bin 
for(i in 2:ncol(ttest.norm.at4)){
temp_data<-ttest.norm.at4[,i]  
temp_out<-t.test(temp_data~ttest.norm.at4$Class)
p.out.at[i]<-temp_out$p.value
}
#create dataframe with rzmt peak names and p-values                           
p.out.at
hist( p.out.at)
peaks2<-colnames(ttest.norm.at4)
View(peaks2)
peaks3<-peaks2[-1]
View(peaks3)
View(p.out.at) 
p.out.at<-p.out.at[-1]                       
p.out.at2<-cbind(peaks3, p.out.at)    
View(p.out.at2)
                           
#### combine ttest p value with feature ranking peaks from SVM RFE ####
#p.out3<-p.out.at2[-1,]                      
#dim(p.out3)
#cbind with ranking list
p.out.at4<-cbind(p.out.at2, featureRankedList[-1])
View(p.out.at4)
class(p.out.at4)
plot(p.out.at4[,3],p.out.at4[,2])
# merge p values with accuracy vector for plot 
#acc1p1<-merge(p.out.at4, acc4merge, by.x="featureRankedList", by.y="no.features")              
#View(acc1p1)

# merge rankings and p values with diffreport for ID
names(XCMSdiff_at_tof_r4_mzrt)                     
p.out.at5<-as.data.frame(p.out.at4)
View(p.out.at5)
SVMrank4ID_at<-merge(p.out.at5, XCMSdiff_at_tof_r4_mzrt, by.x="peaks3", by.y="name")                                            
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/xcms/at_r4_tof_all_fold")
write.csv(SVMrank4ID_at, file="SVMrank4ID_atp071515.csv")  #v3 ranking
                    
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
 
                           
                          