###############
###############
###               CVGLM
###############
###############


### Funkcje składowe:
# 1. Kategoryzująca zadaną kolumną ( i rozszerzenie na ramkę) na podstawie targetu. Zwraca jako rezultat obiekt z inf o tym jak wykonać kategoryzacje na dowolnym zbiorze.
# 2. Funkcja nakładająca na nowy zbiór wyznaczoną kategoryzację.
# 3. Budująca zbior danych do modelowania z CV
# 4. Budująca model i Wyznaczająca zadaną statystykę na zbiorach CV

library(data.table)
library(dplyr)
library(caret)
library(ROCR)

source('CVGLM/ToolFunctionWS.R')



#### Funkcja kategoryzuje i wylicza wartości WOE dla kolumny. Zwraca liste ze slownikami.
categorize.and.woe.feature <- function(x, y, force.factor = FALSE, min.levels.to.continuous = 6 ,eps = 0.1)
{
  ## Test if the x is continuous or factor
  is.contin <- FALSE
  if( force.factor == FALSE && is.numeric(x) && length(unique(x)) >= min.levels.to.continuous ) is.contin <- TRUE
  
  res <- var_partition(variable = x,flag_contin = is.contin, target = y, write_flag = FALSE,eps = eps )
  
  f_var        = res$variable
  t2.counts            = as.table( cbind(tapply(1-y, f_var, FUN = "sum"), tapply(y, f_var, FUN = "sum")) ) 
  colnames(t2.counts)  = c(0,1)
  t2.counts[which(is.na(t2.counts))] = 0                                                                                                      
  t2.prob          = prop.table( t2.counts, margin=2 )
  
  t.bad.prob       = as.table(t2.prob[,2])
  t.good.prob      = as.table(t2.prob[,1])
  
  t.woe            = as.data.frame(log( t.good.prob / t.bad.prob ))
  names(t.woe) <- c('label','woe')
  
  #return( list(cat.dict = res$bands, woe.dict = t.woe))
  colnames( res$bands )[ncol(res$bands)] <- 'label'
  res$bands$label <- as.factor(res$bands$label)
  dplyr::left_join(res$bands, t.woe,by = 'label')
}

set.categorization.and.woe.feature <- function(x, f.trans,flag_label = FALSE)
{
  is.contin <- FALSE
  if( colnames(f.trans)[1] == 'low_cl' ) is.contin <- TRUE
  woe_coding(f.trans ,flag_cont = is.contin, variable = x,flag_label = flag_label)
}


categorize.and.woe.data.frame <- function( data, features, target,min.levels.to.continuous = 6,eps =0.1  )
{
  list.dict <- lapply(features,FUN = function(f){print(f) 
    categorize.and.woe.feature(x = data[,f],y = data[,target],min.levels.to.continuous =min.levels.to.continuous ,eps = eps )} )
  names(list.dict) <- features
  list.dict
}


set.categorization.and.woe.data.frame <- function( data, features, f.trans, flag_label = FALSE)
{
  data.frame(sapply(colnames(data), FUN = function(f){ 
                        if( f %in%features )
                        {
                          set.categorization.and.woe.feature(data[,f],f.trans = f.trans[[f]],flag_label = flag_label)
                        }
                        else{
                          data[,f]
                        }
                }
              )
            )
}


### funkcja ktora zadany zbior dzieli na k-fold, nastepnie 

build.cv.dataset <- function(dataset,features,target,k,seed)
{
  set.seed(seed)
  train.folds <- caret::createFolds(dataset[,target],k = 5,returnTrain = TRUE)
  
  train.cv <- lapply(train.folds, FUN = function(x){
    dataset[x,c(features,target)]
  })
  test.cv <-  lapply(train.folds, FUN = function(x){
    dataset[-x,c(features,target)]
  }) 
  list(train.cv = train.cv,test.cv = test.cv)
}

categorize.and.woe.cv <- function( train.cv,features,target,min.levels.to.continuous = 6,eps =0.1  )
{
  
  dict.cv <-lapply(train.cv, FUN = function(data){
    categorize.and.woe.data.frame(data= data,features = features,target = target, min.levels.to.continuous = min.levels.to.continuous, eps = eps)
  })
  
  
}


set.categorization.and.woe.cv <- function(cv.dataset, cv.trans, features)
{
  k <- length(cv.dataset$train.cv)
  trans.train.cv <- lapply(1:k,FUN = function(i)
    {
    set.categorization.and.woe.data.frame( data = cv.dataset$train.cv[[i]],features = features,f.trans = cv.trans[[i]] )
  })
  trans.test.cv <- lapply(1:k,FUN = function(i)
  {
    set.categorization.and.woe.data.frame( data = cv.dataset$test.cv[[i]],features = features,f.trans = cv.trans[[i]] )
  })
  list(trans.train.cv = trans.train.cv,trans.test.cv=trans.test.cv)
}


build.glm.model.and.predict.cv <- function( trans.cv.dataset, features, target,k)
{
  formula <-as.formula(paste(target,'~',paste(features,collapse = '+')))
  
  formula
  lapply(1:k, FUN = function(i)
         {
    model <- glm(formula,data =  trans.cv.dataset$trans.train.cv[[i]],family = 'binomial')
    predict.trn <- predict(model,trans.cv.dataset$trans.train.cv[[i]])
    predict.tst <- predict(model,trans.cv.dataset$trans.test.cv[[i]])
    list(predict.trn = predict.trn,y.trn =trans.cv.dataset$trans.train.cv[[i]][,target] 
         ,predict.tst=predict.tst,y.tst = trans.cv.dataset$trans.test.cv[[i]][,target])
  })

}


auc.evaluate <- function(predicted,observed)
{
  pred <- ROCR::prediction(predicted,observed)
  perf <- ROCR::performance(pred,"tpr","fpr")
  perf <- ROCR::performance(pred,"auc")
  as.numeric(perf@y.values)
}


auc.evaluate.cv <- function(pred.cv)
{
  sapply(pred.cv, FUN = function(pred){
    auc.evaluate(pred$predict.tst,pred$y.tst)
  })
}




######################
# SAMPLE

data <- iris[,1:4]
data$target <- ifelse(iris$Species =='virginica',1,0)
data$target <- sample(c(1,0),150,r=T)

features <- names(data)[-5]
target <- 'target'
trans.df <- categorize.and.woe.data.frame(data,features,'target',eps = 0.1)

set.categorization.and.woe.data.frame(data = data, features = names(data)[1:2],f.trans = trans.df)


####### TUTAJ MIESO 

cv.dataset <- build.cv.dataset(dataset = data,features = features,target = 'target',k = 5 ,seed = 112)

cv.trans <- categorize.and.woe.cv(train.cv = cv.dataset$train.cv,features = features, target = target)

trans.cv.dataset <- set.categorization.and.woe.cv(cv.dataset,cv.trans,features = features)

pred.cv <- build.glm.model.and.predict.cv(trans.cv.dataset = trans.cv.dataset,features = features[1:2],k = 5,target = 'target')

set.auc.tst <- auc.evaluate.cv(pred.cv)
mean(set.auc.tst)
sd(set.auc.tst)
