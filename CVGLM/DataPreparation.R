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






######################
# SAMPLE

data <- iris[,1:4]
data$target <- ifelse(iris$Species =='virginica',1,0)
data$target <- sample(c(1,0),150,r=T)

trans.df <- categorize.and.woe.data.frame(data,names(data)[-5],'target',eps = 0.00000000001)

set.categorization.and.woe.data.frame(data = data, features = names(data)[1:2],f.trans = trans.df)

