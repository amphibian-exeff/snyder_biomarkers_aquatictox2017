library(e1071)
#### import data processed in MUMA ####
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz2/WorkDir_allmz9000_noutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz2/WorkDir_allmz9000_noutlier/Preprocessing_Data_a/ProcessedTable.csv")
data2<-cbind(class$V1, ProcessedTable)
View(data2)
data2<-data2[-6:-10,]
View(data2)
data2<-data2[-12:-20,]
View(data2)
#data2$Class<-as.factor(data2$class$V1)
data2$Class<-c(1,1,1,1,1,5,5,5,5,5,5)
View(data2)
data2<-data2[,-1:-2]
View(data2)
dim(data2)
#### tuning ####
set.seed(33)
gamma<-c(0.1,1)
tuned <- tune.svm(as.factor(Class)~., data = data2, gamma = gamma, cost = 0.1, 
                  tunecontrol = tune.control(cross = 3), scale=FALSE) 

- sampling method: 3-fold cross validation 
- best parameters:
  gamma cost
0.1  0.1
- best performance: 0.4444444 
- Detailed performance results:
  gamma cost     error dispersion
1   0.1  0.1 0.4444444 0.09622504
2   1.0  0.1 0.4444444 0.09622504

#cost = 10^(-10:10)
#gamma = 0.1

tuned <- tune.svm(as.factor(Class)~., data = data2, gamma = 0.1, cost = 10^(-10:10), 
                  tunecontrol = tune.control(cross = 3), scale=FALSE) 
#was taking a long time

cost.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
gamma.vector<-c(1e-10,1e-5,1e-1,1e1,1e5,1e10)
accuracy.vector<-NULL
m9<-matrix(NA, length(cost.vector), length(gamma.vector))

for (k in seq_along(cost.vector)){
  for (j in seq_along(gamma.vector)){    
    for (i in 1:250){
      model.radial.all <- svm(as.factor(Class)~., data=data2, kernal="radial", gamma=gamma.vector[j], cost=cost.vector[k], scale=FALSE, cross=3)
      accuracy.vector[i]<-model.radial.all$tot.accuracy
    }
    m9[k, j]<-mean(accuracy.vector)  
  }
}
m9
# finished 9 250 SVMs in 9 hours

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

#### run SVM RFE on full feature set ####
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz2/WorkDir_allmz9000_noutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz2/WorkDir_allmz9000_noutlier/Preprocessing_Data_a/ProcessedTable.csv")
X<-ProcessedTable
names(X)
X<-X[-6:-10,]
View(X)
X<-X[-12:-20,]
View(X)
X<-X[,-1]
names(X)
Y<-c(1,1,1,1,1,5,5,5,5,5,5)
Y<-as.factor(Y)
summary(Y)
names(X)
# do RFE
featureRankedList <-svmrfeFeatureRanking(X,Y)
# run SVM with 200 top ranked bins from RFE
svmModel = svm(X[, featureRankedList[1:200]], Y, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
svmModel 
summary(svmModel)
print(svmModel)

#### do 250 bootstrapped leave one out cross validation for SVM with 1-200 bins #### 
accuracy.vector<-NULL
for (i in 1:250){
  svmModel1 = svm(X[, featureRankedList[1:200]], Y, cost = 4, kernel="radial", cross=3 )
  accuracy.vector[i]<-svmModel1$tot.accuracy
}

#### do 250 bootstrapped leave one out cross validation for SVM with 1-9000 bins ####
no.features<-seq(1,9000,by=100)
out.acc<-NULL
for (j in seq_along(no.features)){
  for (i in 1:250){
    svmModel3 = svm(X[, featureRankedList[1:no.features[j]]], Y, cost = 4, kernel="radial", cross=3 )
    accuracy.vector[i]<-svmModel3$tot.accuracy
  }
  out.acc[j]<-mean(accuracy.vector)  
}
out.acc
acc4merge<-cbind(no.features, out.acc)
View(acc4merge)

#### produce figure for paper showing how number of bins included increases classification accuracy ####
setwd("~/Dropbox/amphib_metabolomics/DATA/toad_gt_stage22_round4/results/figures")
png("GT_accuracyperbin.png",height = 4, width = 4, units = 'in', res=600)
no.features<-seq(1,9000,by=100)
plot(out.acc~no.features, las=1, xlab="Number of bins", ylab="Accuracy", ylim=c(0,100), pch=16)
dev.off()

#### zoom in closer to 1 to 250
accuracy.vector2<-NULL
no.features<-seq(1,500,by=20)
out.acc2<-NULL
for (j in seq_along(no.features)){
  for (i in 1:250){
    svmModel3 = svm(X[, featureRankedList[1:no.features[j]]], Y, cost = 4, kernel="radial", cross=3 )
    accuracy.vector2[i]<-svmModel3$tot.accuracy
  }
  out.acc2[j]<-mean(accuracy.vector2)  
}
no.features<-seq(1,500,by=20)
plot(out.acc2~no.features, xlab="Number of features", ylab="Accuracy", ylim=c(0,100), pch=16)
View(out.acc2)

accuracy.vector3<-NULL
no.features<-seq(1,10,by=1)
out.acc3<-NULL
for (j in seq_along(no.features)){
  for (i in 1:250){
    svmModel3 = svm(X[, featureRankedList[1:no.features[j]]], Y, cost = 4, kernel="radial", cross=3 )
    accuracy.vector3[i]<-svmModel3$tot.accuracy
  }
  out.acc3[j]<-mean(accuracy.vector3)  
}
no.features<-seq(1,10,by=1)
plot(out.acc3~no.features, xlab="Number of features", ylab="Accuracy", ylim=c(0,100), pch=16)
View(out.acc3)

#### p value per column  ####
featureRankedList
length(featureRankedList) # 9722
#peaks.rank<-rbind(featureRankedList, X)
#View(peaks.rank)
peaks.t<-t(Y)
View(peaks.t)

# import diffreport
XCMSdiffreport_0416154SVMjoin <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/XCMSdiffreport_0416154SVMjoin.csv")
View(XCMSdiffreport_0416154SVMjoin)
ttest1<-XCMSdiffreport_0416154SVMjoin
dim(ttest1) #9723   48
# delete unnecessary columns
ttest1<-ttest1[,-1]
names(ttest1)
dim(ttest1)
View(ttest1)
ttest1<-ttest1[,-2:-21]
names(ttest1)
# drop 10,50,250 doses
ttest2<-ttest1[,-7:-11]
names(ttest2)
ttest2<-ttest2[,-13:-22]
names(ttest2)
peaks<-ttest2[,1]
ttest3<-ttest2[,-1]
# normalize each column sum column abundance/total
ttest.norm<-matrix(NA, nrow(ttest3), ncol(ttest3))
for (i in 1:ncol(ttest3)){
  outsum<-sum(ttest3[,i])
  for (j in 1:nrow(ttest3)){
    ttest.norm[j,i]<-ttest3[j,i] / outsum
  }
}
View(ttest.norm)
ttest.norm2<-cbind(peaks,as.data.frame(ttest.norm))
View(ttest.norm2)
colnames(ttest.norm2)<-names(ttest2)
View(ttest.norm2)

#transform matrix 
ttest.norm3<-t(ttest.norm2)
View(ttest.norm3)
colnames(ttest.norm3)<-ttest.norm3[1,]
View(ttest.norm3)
ttest.norm3<-ttest.norm3[-1,]
View(ttest.norm3)
Class<-c("c","c","c","c","c","d1250","d1250","d1250","d1250","d1250","d1250")
Class<-as.factor(Class)
class(Class)
ttest.norm4<-cbind(Class, ttest.norm3)
View(ttest.norm4)
ttest.norm4<-as.data.frame(ttest.norm4)
View(ttest.norm4)
class(ttest.norm4)

#convert to numeric
test<-as.numeric(as.character(ttest.norm4[,3]))
class(test)
View(test)
ttest.norm4<-as.numeric((unlist(ttest.norm4[,-1]))                 
for( i in 2:ncol(ttest.norm4)){
  ttest.norm4[,i]<-as.numeric(as.character(ttest.norm4[,i]))
}
summary(ttest.norm4[,5])        
                        
# loop through tttest of each row
t.out.test<-t.test(ttest.norm4[,6]~ttest.norm4$Class )
summary(t.out.test)                        
t.out.test$p.value
#set up empty vector
dim(ttest.norm4)
p.out<-vector()
#loop through each bin 
for(i in 2:ncol(ttest.norm4)){
  temp_data<-ttest.norm4[,i]  
  temp_out<-t.test(temp_data~ttest.norm4$Class)
  p.out[i]<-temp_out$p.value
}
p.out
hist( p.out)
peaks2<-colnames(ttest.norm4)
View(peaks2)
peaks3<-peaks2[-1]
View(peaks3)
View(p.out) 
p.out<-p.out[-1]                       
p.out2<-cbind(peaks3, p.out)    
View(p.out2)

#### combine ttest p value with feature ranking peaks from SVM RFE ####
p.out3<-p.out2[-1,]                      
dim(p.out3)
#cbind with ranking list
p.out4<-cbind(p.out3, featureRankedList)
View(p.out4)
names(p.out4)
class(p.out4)
plot(p.out4[,3],p.out4[,2])
# merge p values with accuracy vector for plot 
acc1p1<-merge(p.out4, acc4merge, by.x="featureRankedList", by.y="no.features")              
View(acc1p1)
# plot
acc1p1[,1]<-as.numeric(as.character(acc1p1[,1]))
acc1p1[,3]<-as.numeric(as.character(acc1p1[,3]))
acc1p1[,4]<-as.numeric(as.character(acc1p1[,4]))                    
View(acc1p1)                                 
acc1p2<-acc1p1[order(acc1p1$featureRankedList),]
View(acc1p2)
names(acc1p2) 
par(mar = c(5, 4, 4, 4) + 0.3)                        
plot(acc1p1$featureRankedList,acc1p1$out.acc, ylim=c(0,100), pch=16, ylab="Accuracy", xlab="Number of features")                     
par(new=TRUE)          
plot(acc1p1$featureRankedList, acc1p1$p.out, axes=FALSE, ann=FALSE, xlab="",ylab="", pch=16, col="red")              
axis(side=4, at = pretty(range(acc1p1$p.out)))
View(acc1p1)
#creat table for ID
# giving merge error
# look for unique?                      
SVMrank4ID<-merge(p.out4, XCMSdiffreport_0416154SVMjoin, by.x="peaks3", by.y="name")                        
View(SVMrank4ID)    
write.csv(SVMrank4ID, file="SVMrank4ID_gtp041615.csv")                        
                        
                        