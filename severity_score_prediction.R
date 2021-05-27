library(ISLR)
library(e1071)
library(caret)

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
merge=data.frame(matrix(0, nrow = 77, ncol = 24657))
row.names(merge)=a$ID
a$train_real_score=a$count*a$scores
View(a)
merge$score=a$train_real_score

for (i in 1:5812){
  for (id in row.names(merge)){
    if (id %in% row.names(train_cd4_gene_expr)){
      merge[which(id==row.names(merge)),i]=train_cd4_gene_expr[which(id==row.names(train_cd4_gene_expr)),i]
    }
    else{
      merge[which(id==row.names(merge)),i]=mean(train_cd4_gene_expr[,i])
    }
  }
}


for (i in 5813:11639){
  for (id in row.names(merge)){
    if (id %in% row.names(train_cd8_gene_expr)){
      merge[which(id==row.names(merge)),i]=train_cd8_gene_expr[which(id==row.names(train_cd8_gene_expr)),i-5812]
    }
    else{
      merge[which(id==row.names(merge)),i]=mean(train_cd8_gene_expr[,i-5812])
    }
  }
}


for (i in 11640:17665){
  for (id in row.names(merge)){
    if (id %in% row.names(train_cd19_gene_expr)){
      merge[which(id==row.names(merge)),i]=train_cd19_gene_expr[which(id==row.names(train_cd19_gene_expr)),i-11639]
    }
    else{
      merge[which(id==row.names(merge)),i]=mean(train_cd19_gene_expr[,i-11639])
    }
  }
}


for (i in 17666:24509){
  for (id in row.names(merge)){
    if (id %in% row.names(train_nasal_gene_expr)){
      merge[which(id==row.names(merge)),i]=train_nasal_gene_expr[which(id==row.names(train_nasal_gene_expr)),i-17665]
    }
    else{
      merge[which(id==row.names(merge)),i]=mean(train_nasal_gene_expr[,i-17665])
    }
  }
}

for (i in 24510:24657){
  for (id in row.names(merge)){
    if (id %in% row.names(train_nasal_microbiome)){
      merge[which(id==row.names(merge)),i]=train_nasal_microbiome[which(id==row.names(train_nasal_microbiome)),i-24509]
    }
    else{
      merge[which(id==row.names(merge)),i]=mean(train_nasal_microbiome[,i-24509])
    }
  }
}

write.csv(merge,"merge.csv")
train.x=data.matrix(merge[,-24568])
train.x=scale(train.x)
train.x[is.nan(train.x)] = 0

# LASSO
set.seed(1)
library(glmnet)
grid=10^seq(-2,2,by = -.1)
lasso.mod=glmnet(data.matrix(merge[,-24658]),merge$score,alpha = 1,lamda=grid)
plot(lasso.mod)
set.seed(1)
cv.out=cv.glmnet(data.matrix(merge[,-24658]),merge$score,alpha = 1)
plot(cv.out)
bestlam=cv.out$lambda.min
lasso.pred=predict(lasso.mod,s=bestlam,newx=data.matrix(merge[,-24658]))
mean((lasso.pred-merge$score)^2)
cv.out$cvm
cv.out$lambda

# prepare test data
test_cd4_gene_expr = as.data.frame(t(read.table("cd4_gene_expr_test.txt")))
test_cd8_gene_expr = as.data.frame(t(read.table("cd8_gene_expr_test.txt")))
test_cd19_gene_expr = as.data.frame(t(read.table("cd19_gene_expr_test.txt")))
test_nasal_gene_expr = as.data.frame(t(read.table("nasal_gene_expr_test.txt")))
test_nasal_microbiome = as.data.frame(t(read.table("nasal_microbiome_test.txt")))

test_ids = sort(Reduce(union,list(row.names(test_cd4_gene_expr),row.names(test_cd8_gene_expr),
                                      row.names(test_cd19_gene_expr),row.names(test_nasal_gene_expr),
                                      row.names(test_nasal_microbiome))))

test=data.frame(matrix(0, nrow = 57, ncol = 24657))
row.names(test)=test_ids

for (i in 1:5812){
  for (id in row.names(test)){
    if (id %in% row.names(test_cd4_gene_expr)){
      test[which(id==row.names(test)),i]=test_cd4_gene_expr[which(id==row.names(test_cd4_gene_expr)),i]
    }
    else{
      test[which(id==row.names(test)),i]=mean(test_cd4_gene_expr[,i])
    }
  }
}


for (i in 5813:11639){
  for (id in row.names(test)){
    if (id %in% row.names(test_cd8_gene_expr)){
      test[which(id==row.names(test)),i]=test_cd8_gene_expr[which(id==row.names(test_cd8_gene_expr)),i-5812]
    }
    else{
      test[which(id==row.names(test)),i]=mean(test_cd8_gene_expr[,i-5812])
    }
  }
}


for (i in 11640:17665){
  for (id in row.names(test)){
    if (id %in% row.names(test_cd19_gene_expr)){
      test[which(id==row.names(test)),i]=test_cd19_gene_expr[which(id==row.names(test_cd19_gene_expr)),i-11639]
    }
    else{
      test[which(id==row.names(test)),i]=mean(test_cd19_gene_expr[,i-11639])
    }
  }
}


for (i in 17666:24509){
  for (id in row.names(test)){
    if (id %in% row.names(test_nasal_gene_expr)){
      test[which(id==row.names(test)),i]=test_nasal_gene_expr[which(id==row.names(test_nasal_gene_expr)),i-17665]
    }
    else{
      test[which(id==row.names(test)),i]=mean(test_nasal_gene_expr[,i-17665])
    }
  }
}

for (i in 24510:24657){
  for (id in row.names(test)){
    if (id %in% row.names(test_nasal_microbiome)){
      test[which(id==row.names(test)),i]=test_nasal_microbiome[which(id==row.names(test_nasal_microbiome)),i-24509]
    }
    else{
      test[which(id==row.names(test)),i]=mean(test_nasal_microbiome[,i-24509])
    }
  }
}
View(test)
write.csv(test,"test.csv")
test.pred=predict(lasso.mod,s=bestlam,newx=data.matrix(test))
testpred=data.frame(test.pred)

prediction_score$severity_score=0
for (id in prediction_score$subject_id){
  if (id %in% gsub("\\X","",row.names(testpred)))
  {
    prediction_score$severity_score[which(id==prediction_score$subject_id)]=
    testpred$X1[which(id==gsub("\\X","",row.names(testpred)))]
  }
  else
  {
    prediction_score$severity_score[which(id==prediction_score$subject_id)]=mean(test.pred)
  }
}

prediction_score$severity_score[prediction_score$severity_score==0]=mean(prediction_score$severity_score[prediction_score$severity_score!=0])
prediction_score$X=NULL
write.csv(prediction_score,"prediction.csv",row.names = FALSE)



