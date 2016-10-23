#############
############# TEST NA OVERFITING
#############


M <- 300
N <- 1000

set.seed(1234)
dataset <- data.frame(matrix(ncol  = M, rnorm(N*M))) 
y <- sample(c(0,1),N,r=T)
dataset <- cbind(dataset, y = y)

features <- paste0('X',1:M)

test.id <- caret::createDataPartition(y,p=0.2,list = F)

train <- dataset[-test.id,]
test <- dataset[test.id,]


######## WERSJA Z FIXED eps

cv.dataset <- build.cv.dataset(train,features = features,target = 'y',k = 10,seed = 12124)

cv.trans <- categorize.and.woe.cv(cv.dataset$train.cv,features,'y')

cv.dataset.trans <- set.categorization.and.woe.cv(cv.dataset,cv.trans, features)



######### zaden model nie zostanie znaleziony - TO jest sukces w tym przypadku
step.auc.cv(cv.dataset.trans,features[1:100],'y',15,10)



####### WERSJA GRID

g.1 <- rep(0.01, length(features))
g.2 <- rep(0.3, length(features))

grid.eps <- data.frame(rbind(g.1,g.2))


res <- stepAUC.cv.grid.eps(train,features = features,target = 'y',k = 10,steps = 5,grid.eps = grid.eps)


tst.auc.mean.eps.grid <- sapply( res, FUN = function(x){ mean(x$resultStepAUC[[2]]$tst.set)})
best.eps.id <- which.max(tst.auc.mean.eps.grid)

( best.eps <- grid.eps[best.eps.id] )
( best.features <- res[[best.eps.id]]$resultStepAUC[[1]] )
( best.auc.tst.set <- res[[best.eps.id]]$resultStepAUC[[2]]$tst.set )
mean(best.auc.tst.set)
sd(best.auc.tst.set)

