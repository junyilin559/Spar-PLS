# Download the Wave Energy Converters dataset from:  
# https://archive.ics.uci.edu/dataset/494/wave+energy+converters  
# Ensure the dataset is stored in the './data' directory.  

Rcpp::sourceCpp("../src/computing.cpp")
source("../src/SparPLS.R", echo=TRUE)
source("../src/RandomPLS.R", echo=TRUE)
source("../src/IFPLS.R", echo=TRUE)
source("../src/MCDPLS.R", echo=TRUE)
source("../src/Utils.R", echo=TRUE)

library(Matrix)
library(MASS)
library(pls)
library(RSpectra)

get_Data = function(loc = 'Perth', num = 49, return_total = FALSE){
  file_name <- paste0('data/WEC_', loc, '_', num, '.csv')
  df <- read.csv(file_name)
  col_std <- apply(df, 2, sd)
  df <- df[, col_std > 0]
  df = scale(df)
  
  p = ncol(df)
  X = df[, 1:(2*num)]
  Y = df[, (2*num+1):p]
  if(return_total){
    Y = Y[, num+2]
  }else{
    Y = Y[, 1:num]
  }
  
  return(list(X = X, Y = Y))
}

data_preprocess = function(Data){
  X = Data$X_train
  p = ncol(X)
  Y = Data$Y_train
  Prob_spar = ElementProb(X, Y)
  Prob_unif = RowProb(X, Y, option = 'unif')
  Prob_lev = RowProb(X, Y, option = 'lev')
  
  Data$Prob_spar = Prob_spar
  Data$Prob_unif = Prob_unif
  Data$Prob_lev = Prob_lev
  return(Data)
}


get_PLS = function(Data, m = 5, r = 500, task = 'Regression'){
  mse = c(0, 0, 0, 0, 0, 0, 0, 0)
  t = c(0, 0, 0, 0, 0, 0, 0, 0)
  
  X_train = Data$X_train
  Y_train = Data$Y_train
  X_test = Data$X_test
  Y_test = Data$Y_test
  
  # PLS
  t1 = Sys.time()
  B <- PLS.coef(X_train, Y_train, m)
  t2 = Sys.time()
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[1] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[1] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  t[1] = as.numeric(difftime(t2, t1, units = "secs"))
  
  # RandomPLS by UNIF
  t1 = Sys.time()
  Prob_unif = RowProb(X_train, Y_train, option = 'unif')
  B <- RandomPLS(X_train, Y_train, m, Prob_unif, r)
  t2 = Sys.time()
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[2] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[2] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  t[2] = as.numeric(difftime(t2, t1, units = "secs"))
  
  # RandomPLS by LEV
  t1 = Sys.time()
  Prob_lev = RowProb(X_train, Y_train, option = 'lev')
  B <- RandomPLS(X_train, Y_train, m, Prob_lev, r)
  t2 = Sys.time()
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[3] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[3] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  t[3] = as.numeric(difftime(t2, t1, units = "secs"))
  
  # IFPLS
  t1 = Sys.time()
  B = IFPLS(X_train, Y_train, m, r)
  t2 = Sys.time()
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[4] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[4] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  t[4] = as.numeric(difftime(t2, t1, units = "secs"))
  
  # MCDPLS
  t1 = Sys.time()
  B = MCDPLS(X_train, Y_train, X_test, m, r)
  t2 = Sys.time()
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[5] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[5] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  t[5] = as.numeric(difftime(t2, t1, units = "secs"))
  
  # SparPLS by xxopt
  t1 = Sys.time()
  Prob_spar = ElementProb(X_train, Y_train)
  B = SparPLS(X_train, Y_train, m, Prob_spar,r)
  t2 = Sys.time()
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[6] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[6] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  t[6] = as.numeric(difftime(t2, t1, units = "secs"))
  
  result = list(mse = mse, t = t)
  return(result)
}


cross_validate = function(Data, m = 5, r = 500, 
                          cv_folds = 10, log = TRUE, 
                          task = 'Regression'){
  mse_fold = matrix(0, cv_folds, 8)
  time_fold = matrix(0, cv_folds, 8)
  
  n = nrow(Data$X)
  folds = sample(rep(1:cv_folds, length.out = n))
  
  for(fold in 1:cv_folds){
    test_idx = which(folds == fold)
    train_idx = setdiff(1:n, test_idx)
    
    X_train = Data$X[train_idx, ]
    Y_train = Data$Y[train_idx, , drop = FALSE]
    X_test = Data$X[test_idx, ]
    Y_test = Data$Y[test_idx, , drop = FALSE]
    
    col_std <- apply(X_train, 2, sd)
    X_train <- X_train[, col_std > 0]
    X_test <- X_test[, col_std > 0]
    
    Data$X_train = X_train
    Data$Y_train = Y_train
    Data$X_test = X_test
    Data$Y_test = Y_test
    Data = data_preprocess(Data)
    
    result = get_PLS(Data, m, r, task)
    
    mse_fold[fold, ] = result$mse
    time_fold[fold, ] = result$t
  }
  
  if(log){
    mse_fold = log(mse_fold)
  }
  
  result_mse = apply(mse_fold, 2, mean)
  result_t = apply(time_fold, 2, mean)
  
  result = list(mse = result_mse, t = result_t)
  return(result)
}

set.seed(1234)
list_results = list()
for(m in c(5, 10, 15, 20)){
  for(loc in c('Perth', 'Sydney')){
    for(num in c(49, 100)){
      Data = get_Data(loc = loc, num = num)
      r = 500
      
      print(paste('loc: ', loc, '; num: ', num, '; m: ', m, '; r: ', r))
      result = cross_validate(Data, m = m, r = r, 
                              cv_folds = 10, log = FALSE, 
                              task = 'Regression')
      print(result)
    }
  }
}

