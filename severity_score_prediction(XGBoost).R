library(ISLR)
library(e1071)
library(glmnet)
library(pls)
library(xgboost)
library(readr)
library(stringr)
library(caret)
library(car)
library(tidyverse) 

train_cd4_gene_expr = as.data.frame(t(read.table("cd4_gene_expr_train.txt")))
train_cd8_gene_expr = as.data.frame(t(read.table("cd8_gene_expr_train.txt")))
train_cd19_gene_expr = as.data.frame(t(read.table("cd19_gene_expr_train.txt")))
train_nasal_gene_expr = as.data.frame(t(read.table("nasal_gene_expr_train.txt")))
train_nasal_microbiome = as.data.frame(t(read.table("nasal_microbiome_train.txt")))
train_severity_score = read.table("severity_score_train.txt",header = TRUE)
prediction_score = read.csv("prediction.csv")

training_ids = sort(Reduce(union,list(row.names(train_cd4_gene_expr),row.names(train_cd8_gene_expr),
  row.names(train_cd19_gene_expr),row.names(train_nasal_gene_expr),
  row.names(train_nasal_microbiome))))
a = NULL
a$ID = training_ids
a=as.data.frame(a)
a$count=0
for (id in a$ID){
  if (id %in% row.names(train_cd4_gene_expr)){
    a$count[which(id == a$ID)]=a$count[which(id == a$ID)]+1
  }
  if (id %in% row.names(train_cd8_gene_expr)){
    a$count[which(id == a$ID)]=a$count[which(id == a$ID)]+1
  }
  if (id %in% row.names(train_cd19_gene_expr)){
    a$count[which(id == a$ID)]=a$count[which(id == a$ID)]+1
  }
  if (id %in% row.names(train_nasal_gene_expr)){
    a$count[which(id == a$ID)]=a$count[which(id == a$ID)]+1
  }
  if (id %in% row.names(train_nasal_microbiome)){
    a$count[which(id == a$ID)]=a$count[which(id == a$ID)]+1
  }
}

for (id in a$ID){
  a$scores[which(id == a$ID)] = 
    train_severity_score$severity_score[which(gsub("\\X","",id)==train_severity_score$subject_id)]/
    a$count[which(id == a$ID)]
}

# assume each of five contributors have same impacts on severity score

train_cd4_gene_expr$score=0
for (id in row.names(train_cd4_gene_expr)){
  train_cd4_gene_expr$score[which(id==row.names(train_cd4_gene_expr))]=a$scores[which(id==a$ID)]
}
train_cd4_gene_expr$score

train_cd8_gene_expr$score=0
for (id in row.names(train_cd8_gene_expr)){
  train_cd8_gene_expr$score[which(id==row.names(train_cd8_gene_expr))]=a$scores[which(id==a$ID)]
}
train_cd8_gene_expr$score

train_cd19_gene_expr$score=0
for (id in row.names(train_cd19_gene_expr)){
  train_cd19_gene_expr$score[which(id==row.names(train_cd19_gene_expr))]=a$scores[which(id==a$ID)]
}
train_cd19_gene_expr$score

train_nasal_gene_expr$score=0
for (id in row.names(train_nasal_gene_expr)){
  train_nasal_gene_expr$score[which(id==row.names(train_nasal_gene_expr))]=a$scores[which(id==a$ID)]
}
train_nasal_gene_expr$score

train_nasal_microbiome$score=0
for (id in row.names(train_nasal_microbiome)){
  train_nasal_microbiome$score[which(id==row.names(train_nasal_microbiome))]=a$scores[which(id==a$ID)]
}
train_nasal_microbiome$score

#SVMs for five different model
# tune.out=tune(svm,score~.,data=train_cd4_gene_expr,kernel="radial",
#              ranges=list(cost=c(.1,1,10,100,1000),
#              gamma=c(.001,.01,.1,.5)))
# summary(tune.out)
svm.best.1=svm(score~.,data=train_cd4_gene_expr,kernel="linear",cost=10)
summary(svm.best.1)
pred.1=predict(svm.best.1,newdata =train_cd4_gene_expr)
mean(pred.1-train_cd4_gene_expr$score)^2

cv_ctrl = trainControl(method = "repeatedcv", repeats = 2,number = 4)
tune_grid = expand.grid(
  nrounds = seq(from = 200, to = 2000, by = 50),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(2, 3, 4),
  gamma = c(5,10,15),
  colsample_bytree = 0.7,
  min_child_weight = 5,
  subsample = 0.8)

xgb_tune = train(score ~.,
                 data=train_cd4_gene_expr,
                 method="xgbTree",
                 metric = "RMSE",
                 objective ='reg:squarederror',
                 trControl=cv_ctrl,
                 tuneGrid=tune_grid)
xgb_tune
xgb_tune$bestTune

dtrain = xgb.DMatrix(data = as.matrix(train_cd4_gene_expr[,-5813]), label=as.matrix(train_cd4_gene_expr$score))
params = list(
  booster = "gbtree",
  objective ='reg:squarederror',
  eta=0.05,
  gamma=5,
  max_depth=2,
  min_child_weight=5,
  subsample=0.8,
  colsample_bytree=0.7
)
as.matrix(train_cd4_gene_expr[,-5813])
library(caret)
test_cd4_gene_expr = as.data.frame(t(read.table("cd4_gene_expr_test.txt")))
pred.1 = predict(xgb_tune$bestTune,newdata=test_cd4_gene_expr)

svm.best.2=svm(score~.,data=train_cd8_gene_expr,kernel="linear",cost=10)
pred.2=predict(svm.best.2,newdata =train_cd8_gene_expr)

svm.best.3=svm(score~.,data=train_cd19_gene_expr,kernel="linear",cost=10)
pred.3=predict(svm.best.3,newdata =train_cd19_gene_expr)

svm.best.4=svm(score~.,data=train_nasal_gene_expr,kernel="linear",cost=10)
pred.4=predict(svm.best.4,newdata=train_nasal_gene_expr)

svm.best.5=svm(score~.,data=train_nasal_microbiome,kernel="linear",cost=10)
pred.5=predict(svm.best.5,newdata=train_nasal_microbiome)

pred1=data.frame(pred.1)
View(pred1)
pred2=data.frame(pred.2)
pred3=data.frame(pred.3)
pred4=data.frame(pred.4)
pred5=data.frame(pred.5)
pred2$pred.2
row.names(pred1)

a$predicted=0
for (id in a$ID){
  if (id %in% row.names(pred1)){
    a$predicted[which(id==a$ID)]=a$predicted[which(id==a$ID)]+
      pred1$pred.1[which(id==row.names(train_cd4_gene_expr))]
  }
  if (id %in% row.names(pred2)){
    a$predicted[which(id==a$ID)]=a$predicted[which(id==a$ID)]+
      pred2$pred.2[which(id==row.names(train_cd8_gene_expr))]
  }
  if (id %in% row.names(pred3)){
    a$predicted[which(id==a$ID)]=a$predicted[which(id==a$ID)]+
      pred3$pred.3[which(id==row.names(train_cd19_gene_expr))]
  }
  if (id %in% row.names(pred4)){
    a$predicted[which(id==a$ID)]=a$predicted[which(id==a$ID)]+
      pred4$pred.4[which(id==row.names(train_nasal_gene_expr))]
  }
  if (id %in% row.names(pred5)){
    a$predicted[which(id==a$ID)]=a$predicted[which(id==a$ID)]+
      pred5$pred.5[which(id==row.names(train_nasal_microbiome))]
  }
}
View(a)

# Test data
test_cd4_gene_expr = as.data.frame(t(read.table("cd4_gene_expr_test.txt")))
View(test_cd4_gene_expr)
test_cd8_gene_expr = as.data.frame(t(read.table("cd8_gene_expr_test.txt")))
View(test_cd8_gene_expr)
test_cd19_gene_expr = as.data.frame(t(read.table("cd19_gene_expr_test.txt")))
View(train_cd19_gene_expr)
test_nasal_gene_expr = as.data.frame(t(read.table("nasal_gene_expr_test.txt")))
View(test_nasal_gene_expr)
test_nasal_microbiome = as.data.frame(t(read.table("nasal_microbiome_test.txt")))
View(test_nasal_microbiome)

test.1=predict(svm.best.1,newdata =test_cd4_gene_expr)
test.2=predict(svm.best.2,newdata =test_cd8_gene_expr)
test.3=predict(svm.best.3,newdata =test_cd19_gene_expr)
test.4=predict(svm.best.4,newdata =test_nasal_gene_expr)
test.5=predict(svm.best.5,newdata =test_nasal_microbiome)

test1=data.frame(test.1)
test2=data.frame(test.2)
test3=data.frame(test.3)
test4=data.frame(test.4)
test5=data.frame(test.5)

gsub("\\X","",row.names(pred1))

prediction_score$severity_score=0
for (id in prediction_score$subject_id){
  if (id %in% gsub("\\X","",row.names(test1))){
    prediction_score$severity_score[which(id==prediction_score$subject_id)]=
      prediction_score$severity_score[which(id==prediction_score$subject_id)]+
      test1$test.1[which(id==gsub("\\X","",row.names(test_cd4_gene_expr)))]
  }
  if (id %in% gsub("\\X","",row.names(test2))){
    prediction_score$severity_score[which(id==prediction_score$subject_id)]=
      prediction_score$severity_score[which(id==prediction_score$subject_id)]+
      test2$test.2[which(id==gsub("\\X","",row.names(test_cd8_gene_expr)))]
  }
  if (id %in% gsub("\\X","",row.names(test3))){
    prediction_score$severity_score[which(id==prediction_score$subject_id)]=
      prediction_score$severity_score[which(id==prediction_score$subject_id)]+
      test3$test.3[which(id==gsub("\\X","",row.names(test_cd19_gene_expr)))]
  }
  if (id %in% gsub("\\X","",row.names(test4))){
    prediction_score$severity_score[which(id==prediction_score$subject_id)]=
      prediction_score$severity_score[which(id==prediction_score$subject_id)]+
      test4$test.4[which(id==gsub("\\X","",row.names(test_nasal_gene_expr)))]
  }
  if (id %in% gsub("\\X","",row.names(test5))){
    prediction_score$severity_score[which(id==prediction_score$subject_id)]=
      prediction_score$severity_score[which(id==prediction_score$subject_id)]+
      test5$test.5[which(id==gsub("\\X","",row.names(test_nasal_microbiome)))]
  }
}
prediction_score$severity_score[prediction_score$severity_score==0]=mean(prediction_score$severity_score[prediction_score$severity_score!=0])
prediction_score$X=NULL
write.csv(prediction_score,"prediction.csv",row.names = FALSE)

mean(a$scores-a$predicted)^2




