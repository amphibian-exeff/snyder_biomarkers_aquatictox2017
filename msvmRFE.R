#### source code for RFE function ####
# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

svmRFE.wrap <- function(test.fold, X, ...) {
  # Wrapper to run svmRFE function while omitting a given test fold
  train.data = X[-test.fold, ]
  test.data  = X[test.fold, ]
  
  # Rank the features
  features.ranked = svmRFE(train.data, ...)
  
  return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
}

svmRFE <- function(X, k=1, halve.above=5000) {
  # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
  n = ncol(X) - 1
  
  # Scale data up front so it doesn't have to be redone each pass
  cat('Scaling data...')
  X[, -1] = scale(X[, -1])
  cat('Done!\n')
  flush.console()
  
  pb = txtProgressBar(1, n, 1, style=3)
  
  i.surviving = 1:n
  i.ranked    = n
  ranked.list = vector(length=n)
  
  # Recurse through all the features
  while(length(i.surviving) > 0) {
    if(k > 1) {
      # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
      folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
      folds = lapply(1:k, function(x) which(folds == x))
      
      # Obtain weights for each training set
      w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
      w = do.call(rbind, w)
      
      # Normalize each weights vector
      w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
      
      # Compute ranking criteria
      v    = w * w
      vbar = apply(v, 2, mean)
      vsd  = apply(v, 2, sd)
      c    = vbar / vsd
    } else {
      # Only do 1 pass (i.e. regular SVM-RFE)
      w = getWeights(NULL, X[, c(1, 1+i.surviving)])
      c = w * w
    }
    
    # Rank the features
    ranking = sort(c, index.return=T)$ix
    if(length(i.surviving) == 1) {
      ranking = 1
    }
    
    if(length(i.surviving) > halve.above) {
      # Cut features in half until less than halve.above
      nfeat = length(i.surviving)
      ncut  = round(nfeat / 2)
      n     = nfeat - ncut
      
      cat('Features halved from', nfeat, 'to', n, '\n')
      flush.console()
      
      pb = txtProgressBar(1, n, 1, style=3)
      
    } else ncut = 1
    
    # Update feature list
    ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
    i.ranked    = i.ranked - ncut
    i.surviving = i.surviving[-ranking[1:ncut]]
    
    setTxtProgressBar(pb, n-length(i.surviving))
    flush.console()
  }
  
  close(pb)
  
  return (ranked.list)
}

getWeights <- function(test.fold, X) {
  # Fit a linear SVM model and obtain feature weights
  train.data = X
  if(!is.null(test.fold)) train.data = X[-test.fold, ]
  
  svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
                 scale=F, type="C-classification", kernel="radial")
  
  t(svmModel$coefs) %*% svmModel$SV
}

WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
  # Compile feature rankings across multiple folds
  featureID = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$ix
  avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
  feature.name = colnames(input[, -1])[featureID]
  features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
  if(save==T) {
    write.table(features.ranked, file=file, quote=F, row.names=F)
  } else {
    features.ranked
  }
}

FeatSweep.wrap <- function(i, results, input) {
  # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
  svm.list = lapply(results, function(x) tune(svm,
                                              train.x      = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                              train.y      = input[x$train.data.ids, 1],
                                              validation.x = input[x$test.data.ids, 1+x$feature.ids[1:i]],
                                              validation.y = input[x$test.data.ids, 1],
                                              # Optimize SVM hyperparamters
                                              ranges       = tune(svm,
                                                                  train.x = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                                                  train.y = input[x$train.data.ids, 1],
                                                                  ranges  = list(gamma=2^(-12:0), cost=2^(-6:6)))$best.par,
                                              tunecontrol  = tune.control(sampling='fix'))$perf)
  
  error = mean(sapply(svm.list, function(x) x$error))
  return(list(svm.list=svm.list, error=error))
}

PlotErrors <- function(errors, errors2=NULL, no.info=0.5, ylim=range(c(errors, errors2), na.rm=T), xlab='Number of Features',  ylab='10x CV Error') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='black') {
    lines(which(!is.na(errors)), na.omit(x), col=col)
    points(which.min(x), min(x, na.rm=T), col='red')
    text(which.min(x), min(x, na.rm=T), paste(which.min(x), '-', format(min(x, na.rm=T), dig=3)), pos=4, col='red', cex=0.75)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=3)
}

#### model test c vs. 1250 for subset of m/z values ####
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_Subset/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_Subset/Preprocessing_Data_a/ProcessedTable.csv")
data2<-cbind(class$V1, ProcessedTable)
View(data2)
data2<-data2[-6:-9,]
data2<-data2[-16:-20,]
data2<-data2[-12:-15,]
View(data2)
data2<-data2[,-2]
#data2$Class<-as.factor(data2$class$V1)
data2$Class<-c(1,1,1,1,1,5,5,5,5,5,5)
data2<-data2[,-1]
data2$Class<-as.factor(data2$Class)
data2<-data2[,-131]
names(data2)
class<-data2$Class
data3<-cbind(class, data2)
data3<-data3[,-131]
names(data3)
test1<-svmRFE(data3)
test1
[1]  77  13  44  50  59 125 101   2  47  66  73  61   4 111  20  45  93 128
[19]  85  46  91  43 127 123  69  72 114  12  65  98  17  42  78 118  39  71
[37]  87  19  80  53   3  15 117  81  55  22  35  84  27  70 113 100  49  68
[55]  97  79 112  86  56  76  60  41 108 122  88  74   6  63  96  32  25  29
[73]  58  37   8  92 102  67  36  38 115  64  40 120 104   1  48  99 129  26
[91] 110 109  83  90  31  82  30   9  62  16   7  75  23  89 103  33  28  54
[109] 105  57  24   5  11 121  21  18 116 126  94  95  51 107  52  10 124 119
[127]  34 106  14
 

model <- svm(data3[,test1[1:40]],  kernal="linear", cross=3, scale=FALSE)
print(model)
summary(model)

#### test rfe svm for all m/z values #### 
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz_noutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz_noutlier/Preprocessing_Data_a/ProcessedTable.csv")
data2<-cbind(class$V1, ProcessedTable)
View(data2)
data2<-data2[-6:-10,]
data2<-data2[-12:-20,]
View(data2)
data2<-data2[,-2]
#data2$Class<-as.factor(data2$class$V1)
data2$Class<-c(1,1,1,1,1,5,5,5,5,5,5)
class<-c(1,1,1,1,1,5,5,5,5,5,5)
names(data2)
data2$Class<-as.factor(data2$Class)
data2<-data2[,-1010]
names(data2)
class<-data2$Class
data3<-cbind(class, data2)
data3<-data3[,-1011]
names(data3)
data3<-data3[,-2]
View(data3)
test1<-svmRFE(data3,k=1)
test1
model <- svm(data3[,test1[1:20]],  kernal="radial", cross=4, scale=FALSE)
print(model)
summary(model)

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

#### test of Guyon implementation #### 
class <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz_noutlier/Preprocessing_Data_a/class.csv")
ProcessedTable <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/WorkDir_allmz_noutlier/Preprocessing_Data_a/ProcessedTable.csv")
Y<-ProcessedTable
names(Y)
Y<-Y[-6:-10,]
View(Y)
Y<-Y[-12:-20,]
View(Y)
Y<-Y[,-1]
names(Y)
#data2$Class<-as.factor(data2$class$V1)
X<-c(1,1,1,1,1,5,5,5,5,5,5)
X<-as.factor(X)
summary(X)
featureRankedList <-svmrfeFeatureRanking(Y,X)
svmModel = svm(Y[, featureRankedList[1:200]], X, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
svmModel 
summary(svmModel)
print(svmModel)

accuracy.vector<-NULL
for (i in 1:250){
svmModel1 = svm(Y[, featureRankedList[1:200]], X, cost = 4, kernel="radial", cross=3 )
accuracy.vector[i]<-svmModel1$tot.accuracy
}
mean(accuracy.vector)

#take  5% , 10% 20% , to 100% and then do SVM w/ cross fold validation to get error rate
#no.features<-c(10,100,200,300,400,500,600,700,800,900,1000)
no.features<-c(1,10,100,200,300,400,500,600,700,800,900,1000)
out.acc<-NULL
for (j in seq_along(no.features)){
  for (i in 1:250){
    svmModel1 = svm(Y[, featureRankedList[1:no.features[j]]], X, cost = 4, kernel="radial", cross=3 )
    accuracy.vector[i]<-svmModel1$tot.accuracy
}
 out.acc[j]<-mean(accuracy.vector)  
}
out.acc
plot(out.acc~no.features)
class(Y)
colnames(Y[,874])

#### try again without #1 ranked column ####
Y2<-Y[,-874]
featureRankedList <-svmrfeFeatureRanking(Y2,X)
svmModel = svm(Y2[, featureRankedList[1:200]], X, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
svmModel 
summary(svmModel)
print(svmModel)

#no.features<-c(1,10,100,200,300,400,500,600,700,800,900,1000)
no.features<-c(1,10,100,200)
out.acc2<-NULL
accuracy.vector2<-NULL
for (j in seq_along(no.features)){
  for (i in 1:250){
    svmModel2 = svm(Y2[, featureRankedList[1:no.features[j]]], X, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
    accuracy.vector2[i]<-svmModel2$tot.accuracy
  }
  out.acc2[j]<-mean(accuracy.vector2)  
}
out.acc2
plot(out.acc2~no.features)

#### 0 to 1000 features by 5 step ####
featureRankedList <-svmrfeFeatureRanking(Y,X)
no.features<-seq(1,1000,by=10)
out.acc3<-NULL
accuracy.vector3<-NULL
for (j in seq_along(no.features)){
  for (i in 1:250){
    svmModel3 = svm(Y[, featureRankedList[1:no.features[j]]], X, cost = 6.4e9, gamma=1.0e-10, kernel="radial", cross=3 )
    accuracy.vector3[i]<-svmModel3$tot.accuracy
  }
  out.acc3[j]<-mean(accuracy.vector3)  
}
out.acc3
plot(out.acc3~no.features, xlab="Number of Features", ylab="Accuracy", pch=16, ylim=c(0,100))

#### p value per column plot by rank SVM score ####
featureRankedList
peaks.rank<-rbind(featureRankedList, Y)
View(peaks.rank)
peaks.t<-t(Y)
View(peaks.t)

# import diffreport
XCMSdiffreport_0416154SVMjoin <- read.csv("~/Dropbox/amphib_metabolomics/DATA/gt_stage22_TOF041515/xcms_polar/results/XCMSdiffreport_0416154SVMjoin.csv")
View(XCMSdiffreport_0416154SVMjoin)
ttest1<-XCMSdiffreport_0416154SVMjoin
# delete unnecessary columns
ttest1<-ttest1[,-1]
names(ttest1)
ttest1<-ttest1[,-2:-21]
names(ttest1)
# drop 10,50,250 doses
ttest2<-ttest1[,-7:-11]
names(ttest2)
ttest2<-ttest2[,-13:-22]
names(ttest2)
# normalize each column sum column abundance/total
ttest.norm<-matrix(NA, nrow(ttest2), ncol(ttest2))
#View(ttest.norm)
for (i in 2:ncol(ttest2)){
  outsum<-sum(ttest2[,i])
  for (j in 1:nrow(ttest2)){
    ttest.norm[j,i]<-ttest2[j,i] / outsum
  }
}
View(ttest.norm)
# loop through tttest of each row

# combine feature ranking with peaks

# transpose to long format

# merge with ranking list


#### cross species comparison of metabolomic profiles ####
#### ROC R package ROCR #### 





