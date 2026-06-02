Rcpp::sourceCpp("../src/computing.cpp")
source("../src/SparPLS.R", echo=TRUE)
source("../src/RandomPLS.R", echo=TRUE)

library(Matrix)
library(MASS)
library(pls)
library(RSpectra)

random_orthogonal_matrix <- function(dim) {
  Q <- qr.Q(qr(matrix(rnorm(dim * dim), nrow = dim)))
  return(Q)
}

Generate_Data = function(n = 50, p = 500, l = 100, m = 5, R = 0.5){
  P <- random_orthogonal_matrix(p)[1:m, ]
  Q <- random_orthogonal_matrix(l)[1:m, ]
  beta <- t(P) %*% solve(P %*% t(P)) %*% Q
  Sigma <- matrix(0, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      Sigma[i, j] <- 0.5^(abs(i - j))
    }
  }
  eigen_decomp <- eigen(Sigma)
  eigenvalues <- eigen_decomp$values
  eigenvectors <- eigen_decomp$vectors
  Sigma_1_2 <- eigenvectors %*% diag(sqrt(eigenvalues)) %*% t(eigenvectors)
  
  Scale <- diag(exp((0:(m-1)) * log(20) / (m-1) / 2))
  SNR_X <- R
  SNR_Y <- R
  
  T_mat <- mvrnorm(n = n, mu = rep(0, m), Sigma = diag(m))
  X0 <- T_mat %*% P
  Y0 <- T_mat %*% Q
  
  X_std <- sd(X0)
  Y_std <- sd(Y0)
  
  E_mat <- mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p) * (X_std / sqrt(SNR_X))^2)
  F_mat <- mvrnorm(n = n, mu = rep(0, l), Sigma = diag(l) * (Y_std / sqrt(SNR_Y))^2)
  
  X <- T_mat %*% P + E_mat
  Y <- T_mat %*% Q + F_mat
  return(list(X = X, Y = Y, beta = beta, P = P, Q = Q))
}

data_preprocess = function(Data){
  X = Data$X_train
  p = ncol(X)
  Y = Data$Y_train
  Prob_spar = ElementProb(X, Y)
  Prob_unif = RowProb(X, Y, option = 'unif')
  Prob_lev = RowProb(X, Y, option = 'lev')
  
  Data$Prob_spar = Prob_spar
  Data$Prob_lev = Prob_lev
  Data$Prob_unif = Prob_unif
  return(Data)
}

get_PLS = function(Data, m = 5){
  bmse = c(0, 0, 0, 0)
  qmse = c(0, 0, 0, 0)
  pmse = c(0, 0, 0, 0)
  
  X_train = Data$X_train
  Y_train = Data$Y_train
  X_test = Data$X_test
  Y_test = Data$Y_test
  
  n = nrow(X_train)
  r = round(n * 0.1)
  
  # PLS
  B <- PLS.coef(X_train, Y_train, m)
  Y_pred = X_test %*% B
  Q_diff = Data$P %*% (B - Data$beta)
  bmse[1] = mean((B - Data$beta)*(B - Data$beta))
  qmse[1] = mean(Q_diff*Q_diff)
  pmse[1] = mean((Y_pred - Y_test)^2) 
  
  # RandomPLS by UNIF
  B <- RandomPLS(X_train, Y_train, m, Data$Prob_unif, r*4)
  Y_pred = X_test %*% B
  Q_diff = Data$P %*% (B - Data$beta)
  bmse[2] = mean((B - Data$beta)*(B - Data$beta))
  qmse[2] = mean(Q_diff*Q_diff)
  pmse[2] = mean((Y_pred - Y_test)^2) 
  
  # RandomPLS by LEV
  B <- RandomPLS(X_train, Y_train, m, Data$Prob_lev, r*4)
  Y_pred = X_test %*% B
  Q_diff = Data$P %*% (B - Data$beta)
  bmse[3] = mean((B - Data$beta)*(B - Data$beta))
  qmse[3] = mean(Q_diff*Q_diff)
  pmse[3] = mean((Y_pred - Y_test)^2) 
  
  # SparPLS-0.1
  B = SparPLS(X_train, Y_train, m, Data$Prob_spar, r)
  Y_pred = X_test %*% B
  Q_diff = Data$P %*% (B - Data$beta)
  bmse[4] = mean((B - Data$beta)*(B - Data$beta))
  qmse[4] = mean(Q_diff*Q_diff)
  pmse[4] = mean((Y_pred - Y_test)^2) 
  
  # SparPLS-0.2
  B = SparPLS(X_train, Y_train, m, Data$Prob_spar, r*2)
  Y_pred = X_test %*% B
  Q_diff = Data$P %*% (B - Data$beta)
  bmse[5] = mean((B - Data$beta)*(B - Data$beta))
  qmse[5] = mean(Q_diff*Q_diff)
  pmse[5] = mean((Y_pred - Y_test)^2) 
  
  # SparPLS-0.4
  B = SparPLS(X_train, Y_train, m, Data$Prob_spar, r*4)
  Y_pred = X_test %*% B
  Q_diff = Data$P %*% (B - Data$beta)
  bmse[6] = mean((B - Data$beta)*(B - Data$beta))
  qmse[6] = mean(Q_diff*Q_diff)
  pmse[6] = mean((Y_pred - Y_test)^2) 
  
  
  return(list(bmse = bmse, qmse = qmse, pmse = pmse))
}

cross_validate = function(Data, m = 5, 
                          cv_folds = 10, log_out = TRUE){
  bmse_fold = matrix(0, cv_folds, 6)
  qmse_fold = matrix(0, cv_folds, 6)
  pmse_fold = matrix(0, cv_folds, 6)
  
  n = nrow(Data$X)
  folds = sample(rep(1:cv_folds, length.out = n))
  
  for(fold in 1:cv_folds){
    test_idx = which(folds == fold)
    train_idx = setdiff(1:n, test_idx)
    
    X_train = Data$X[train_idx, ]
    Y_train = Data$Y[train_idx, ]
    X_test = Data$X[test_idx, ]
    Y_test = Data$Y[test_idx, ]
    
    Data$X_train = X_train
    Data$Y_train = Y_train
    Data$X_test = X_test
    Data$Y_test = Y_test
    Data = data_preprocess(Data)
    
    result = get_PLS(Data, m)
    
    bmse_fold[fold, ] = result$bmse
    qmse_fold[fold, ] = result$qmse
    pmse_fold[fold, ] = result$pmse
  }
  
  if(log_out){
    bmse_fold = log(bmse_fold)
    qmse_fold = log(qmse_fold)
    pmse_fold = log(pmse_fold)
  }
  
  result_bmse = apply(bmse_fold, 2, mean)
  result_qmse = apply(qmse_fold, 2, mean)
  result_pmse = apply(pmse_fold, 2, mean)
  
  result = list(result_bmse = result_bmse, result_qmse = result_qmse, result_pmse = result_pmse)
  return(result)
}

set.seed(1234)
for(m in c(5, 10, 15)){
  Data = Generate_Data(m = m)
  result = cross_validate(Data, m = m, cv_folds = 10, log = TRUE)
  print(result)
}
