##########
########## Przygotowujemy dane testowe: PRzypadek 5 cech zależnych
#########
M <- 300
N <- 1000

set.seed(1234)
dataset <- data.frame(matrix(ncol  = M, rnorm(N*M))) 
real.feature <- sample(paste0('X',1:30), 5)
real.feature
real.coef <- rnorm(mean = 0, sd = 2, 5)
real.coef

y.surogate <- rowSums(t( t(dataset[,real.feature])*real.coef))
y.surogate
sd(y.surogate)
y.surogate.rnd <- y.surogate + rnorm(N,mean = 0, sd = 1)


y <- y.surogate.rnd < mean(y.surogate.rnd)

dataset <- cbind(dataset, y = y)

features <- paste0('X',1:M)
#########################

# wybieramy zbiór testowy

test.id <- caret::createDataPartition(y,p=0.2,list = F)

train <- dataset[-test.id,]
test <- dataset[test.id,]


######################
###### Z ciekawosci zobaczmy model jaki by powstał tylko na skorelowanych zmiennych
(real.formula <- paste('y ~0+',paste(real.feature,collapse = ' + ')))

model.real <- glm(real.formula, data = train, family = 'binomial')

coef(model.real)
real.coef


trn.pred <- predict(model.real,train,type = 'response' )
tst.pred <- predict(model.real,test, type = 'response' )

### skutecznosc z CUT OFF 0.5 ( zbior zbalansowany)
tst.shot <- tst.pred > 0.5
t <- table(tst.shot, test$y)
rownames(t) <- c(0,1)
colnames(t) <- c(0,1)
confusionMatrix(t)

auc.evaluate(trn.pred,train$y)
auc.evaluate(tst.pred,test$y)



############  TERAZ SPRAWDZMY JAKIE  MIAŁBY statystyki
############  model zlozony z korelowanych z targetem cech wyznacozny na CV Folds
########
######## Zakładamy dość dokladną kategoryzację z eps = 0.01
cv.dataset <- build.cv.dataset(train,features = features,target = 'y',k = 10,seed = 12124)

cv.trans <- categorize.and.woe.cv(cv.dataset$train.cv,features,'y',eps = 0.01)

cv.dataset.trans <- set.categorization.and.woe.cv(cv.dataset,cv.trans, features)


###### REAL 
pred.cv <- build.glm.model.and.predict.cv(cv.dataset.trans,real.feature,target = 'y',k = 10)

auc.set <- auc.evaluate.cv(pred.cv)
auc.set
mean(auc.set$tst.set)
sd(auc.set$tst.set)

mean(auc.set$trn.set)
sd(auc.set$trn.set)

######  Warto to porownac ze statystykami na zbiorze testowym wczesniej wyliczonymi

##############################

####  A Teraz uruchamiamy armatę 
#####
######### GRID Search 

######  Wezmy tylko dwie kategoryzacje i wszystkie cechy tak samo kategorzyzowane

#### Korzystamy tu wprost ze zbioru train  Wszystko sie robi \automatycznie
g.1 <- rep(0.01, length(features))
g.2 <- rep(0.3, length(features))

grid.eps <- data.frame(rbind(g.1,g.2))


res <- stepAUC.cv.grid.eps(train,features = features,target = 'y',k = 10,steps = 5,grid.eps = grid.eps)


tst.auc.mean.eps.grid <- sapply( res, FUN = function(x){ mean(x$resultStepAUC[[2]]$tst.set)})

( best.eps.id <- which.max(tst.auc.mean.eps.grid) )
( best.eps <- grid.eps[best.eps.id] )
( best.features <- res[[best.eps.id]]$resultStepAUC[[1]] )
( best.auc.tst.set <- res[[best.eps.id]]$resultStepAUC[[2]]$tst.set )
mean(best.auc.tst.set)
sd(best.auc.tst.set)


