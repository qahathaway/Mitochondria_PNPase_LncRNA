library(tidyverse)
library(reshape2)
library(xgboost)
library(randomForest)
library(rfUtilities)

Features_All<-read.csv("path_to_your_files.csv", row.names=1)
Features_All[,2] <- as.numeric(as.character(Features_All[,2]))
Features_All[,1] <- as.factor(as.character(Features_All[,1]))

head(Features_All)
summary(Features_All)



###Classification Trees###
# Classification CV
Features_All$label <- as.factor(Features_All$label)
set.seed(100)
( rf.mdl <- randomForest(Features_All[,2:135],Features_All[,"label"], ntree=5000) )
( rf.cv <- rf.crossValidation(rf.mdl, Features_All[,2:135], p=0.10, n=10, ntree=5000) )


# Plot cross validation verses model producers accuracy
par(mfrow=c(2,2))
plot(rf.cv, type = "cv", main = "CV producers accuracy")
plot(rf.cv, type = "model", main = "Model producers accuracy")

# Plot cross validation verses model oob
plot(rf.cv, type = "cv", stat = "oob", main = "CV oob error")
plot(rf.cv, type = "model", stat = "oob", main = "Model oob error")


# Classification Single

set.seed(100)
ind = sample(2, nrow(Features_All), replace=TRUE, prob=c(0.7,0.3))
trainData = Features_All_kmer50.1[ind==1,]
testData = Features_All_kmer50.1[ind==2,]

Features_All_rf = randomForest(label~., data=trainData, ntree=5000, proximity=TRUE)
table(predict(Features_All_rf), trainData$label)
Features_All_rf
par(mfrow=c(1,1))
plot(Features_All_rf)
importance(Features_All_rf)

Features_All_Pred = predict(Features_All_rf, newdata=testData)
table(Features_All_Pred, testData$label)

plot(margin(Features_All_rf, testData$label))
CM = table(Features_All_Pred, testData$label)

accuracy = (sum(diag(CM)))/sum(CM)
accuracy

library(ROCR)
predictions=as.vector(Features_All_rf$votes[,2])
pred=prediction(predictions,trainData$label)

perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]

perf_ROC=performance(pred,"tpr","fpr")
plot(perf_ROC, main="ROC plot")
text(0.5,0.3,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
