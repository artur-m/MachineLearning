##########
########## where are predictors :)
##########



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
  trn.set <- sapply(pred.cv, FUN = function(pred){
    auc.evaluate(pred$predict.trn,pred$y.trn)
  })
  
  tst.set <- sapply(pred.cv, FUN = function(pred){
    auc.evaluate(pred$predict.tst,pred$y.tst)
  })
  list(trn.set = trn.set, tst.set = tst.set)
}


step.auc.cv <- function( cv.dataset, features, target, steps, k)
{
  selected.features <- c()
  current.selected.features <- c()
  avaiable.features <- features
  formula <- as.formula(paste(target,'~.'))
  best.auc  <- random.auc.set.cv(k)
  best.s.auc <- random.auc.set.cv(k)
  s <- 1
  while( s < steps)
  {
    best.s.auc <- random.auc.set.cv(k)
    best.s.feature <- -1
    for( f in avaiable.features)
    {
      current.selected.features <- c(selected.features,f)
      pred.cv <- build.glm.model.and.predict.cv(trans.cv.dataset = cv.dataset
                                                ,features = current.selected.features
                                                ,target = target,k = k)
      
      current.auc.set <- auc.evaluate.cv(pred.cv)
      if( test.predict.ability.cv( current.auc.set) )
      {
        if( auc.model.cmp.cv(best.s.auc, current.auc.set)  )
        {
          best.s.auc = current.auc.set
          best.s.feature <- f
        }
        
      }
      
      
    }
    print(best.s.feature)
    if( auc.model.cmp.cv(best.auc, best.s.auc))
    {
      avaiable.features <- avaiable.features[best.s.feature!=avaiable.features]
      selected.features <- c(selected.features,best.s.feature)
      s <- s+1
      best.auc <- best.s.auc
      print(selected.features)
      print(best.s.auc)
      
    }
    else
    {
      break;
    }
    
    
    
  }
  list(selected.features, best.auc)
}

test.predict.ability.cv <- function( set.auc.cv)
{
  if( min( abs(set.auc.cv$trn.set - set.auc.cv$tst.set) ) > 0.1 ) return(FALSE)
  sd.tst <- sd(set.auc.cv$tst)
  sd.trn <- sd(set.auc.cv$trn)
  if( sd.tst > 0.1 || sd.trn > 0.1 ) return(FALSE)
  if( any( set.auc.cv$tst.set - sd.tst < 0.5 ) || any( set.auc.cv$trn.set - sd.trn < 0.5 ) ) return(FALSE)
  return(TRUE)
}

### x < y
auc.model.cmp.cv <- function( set.auc.cv.x, set.auc.cv.y)
{
  mean(set.auc.cv.x$tst.set) < mean(set.auc.cv.y$tst.set)
}

random.auc.set.cv <- function(k)
{
  list( trn.set = rep(0.5,k),tst.set = rep(0.5,k))
}


### petla po eps w kategoryzacji 

stepAUC.cv.grid.eps <- function( train, features, target, steps, k, grid.eps,seed = 12124)
{
  
  cv.dataset <- build.cv.dataset(train,features = features,target = target,k = k,seed = seed)
  result.stepAUC.list <- list()
  for( i in 1:nrow(grid.eps))
  {
    cv.trans <- categorize.and.woe.cv(cv.dataset$train.cv,features,'y',eps = as.numeric(grid.eps[i,]))
    cv.dataset.trans <- set.categorization.and.woe.cv(cv.dataset,cv.trans, features)
    
    resStepAUC <- step.auc.cv(cv.dataset.trans,features,target,steps = steps,k)
    result.stepAUC.list[[i]] <- list()
    result.stepAUC.list[[i]]$trans <- cv.trans 
    result.stepAUC.list[[i]]$resultStepAUC <- resStepAUC
    
  }
  
  return(result.stepAUC.list)
}

##########