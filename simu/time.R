Rcpp::sourceCpp("../src/computing.cpp")
source("../src/SparPLS.R", echo=TRUE)
source("../src/RandomPLS.R", echo=TRUE)
source("../src/IFPLS.R", echo=TRUE)
source("../src/MCDPLS.R", echo=TRUE)

library(Matrix)
library(MASS)
library(pls)
library(RSpectra)

random_orthogonal_matrix <- function(dim) {
  Q <- qr.Q(qr(matrix(rnorm(dim * dim), nrow = dim)))
  return(Q)
}

Generate_Data = function(n = 50000, p = 100, l = 50, m = 5, R = 0.25){
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
  diag_XX = diag(t(X) %*% X)
  return(list(X = X, Y = Y, diag_XX = diag_XX, beta = beta, P = P, Q = Q))
}

data_preprocess = function(Data){
  X = Data$X
  p = ncol(X)
  Y = Data$Y
  Prob_spar = ElementProb(X, Y)
  Prob_unif = RowProb(X, Y, option = 'unif')
  Prob_lev = ElementProb(X, Y)
  Prob_lev = Prob_lev * (Data$diag_XX/min(Data$diag_XX))
  Prob_lev = matrix(1, nrow(X), p)
  Prob_lev = Prob_lev / sum(Prob_lev)
  
  Data$Prob_spar = Prob_spar
  Data$Prob_lev = Prob_lev
  Data$Prob_unif = Prob_unif
  return(Data)
}

get_PLS = function(Data, m = 5, r = 500){
  t = c(0,0,0,0)
  
  # PLS
  t1 = Sys.time()
  B <- PLS.coef(Data$X, Data$Y, m)
  t2 = Sys.time()
  t[1] = t[1] + as.numeric(difftime(t2, t1, units = "secs"))
  
  # RandomPLS by UNIF
  t1 = Sys.time()
  Prob_unif = RowProb(Data$X, Data$Y, option = 'unif')
  B <- RandomPLS(Data$X, Data$Y, m, Data$Prob_unif, r)
  t2 = Sys.time()
  t[2] = t[2] + as.numeric(difftime(t2, t1, units = "secs"))
  
  # RandomPLS by LEV
  t1 = Sys.time()
  Prob_lev = RowProb(Data$X, Data$Y, option = 'lev')
  B <- SparPLS(Data$X, Data$Y, m, Data$Prob_lev, r)
  t2 = Sys.time()
  t[3] = t[3] + as.numeric(difftime(t2, t1, units = "secs"))
  
  # SparPLS by xxopt
  t1 = Sys.time()
  n0 = nrow(Data$X)
  p0 = ncol(Data$X)
  Prob_spar  = matrix(1/n0/p0, n0,p0)
  B = SparPLS(Data$X, Data$Y, m, Data$Prob_spar, r)
  t2 = Sys.time()
  t[4] = t[4] + as.numeric(difftime(t2, t1, units = "secs"))
  
  return(t)
}

get_PLS_add = function(Data, m = 5, r = 500){
  t = c(0,0)
  # IFPLS
  t1 = Sys.time()
  B = IFPLS(Data$X, Data$Y, m, r)
  t2 = Sys.time()
  t[1] = t[1] + as.numeric(difftime(t2, t1, units = "secs"))
  
  # MCDPLS
  t1 = Sys.time()
  B = MCDPLS(Data$X, Data$Y, Data$X, m, r, simu = TRUE)
  t2 = Sys.time()
  t[2] = t[2] + as.numeric(difftime(t2, t1, units = "secs"))
  
  return(t)
}

repeat_experiments = function(Data, m = 5, r = 500, 
                              rep = 10, log = FALSE){
  result_mat = matrix(0, rep, 6)
  for(i in c(1:rep)){
    result1 = get_PLS(Data, m, r)
    result2 = get_PLS_add(Data, m, r)
    result_mat[i, 1:4] = result1
    result_mat[i, 5:6] = result2
  }
  if(log){
    print(result_mat)
    result_mat = log(result_mat)
    print(result_mat)
  }
  result_mean = apply(result_mat, 2, mean)
  result_std = apply(result_mat, 2, sd)
  result = list(mean = result_mean, 
                std = result_std)
  return(result)
}

df_init = function(subsample_sizes = c(500, 1000, 1500, 2000, 2500),
                   Method = c("FULL", "UNIF", "LEV", "SPAR", "IF", "MCD"),
                   M_list = c("m = 5", "m = 10")){
  set.seed(123)
  data <- expand.grid(
    Method = Method,
    subsample = subsample_sizes,
    M = M_list
  )
  
  data$MSE <- rnorm(nrow(data), mean = 10, sd = 2)
  data$MSE_error <- rnorm(nrow(data), mean = 1, sd = 0.5)
  data$subplot_label <- paste(data$R, data$M, sep = ",")
  return(data)
}



set.seed(1234)
df = df_init()
k0 = 1
R = 1
for(m in c(5, 10)){
  Data = Generate_Data(m = m, R = R)
  Data = data_preprocess(Data)
  for(r in c(500, 1000, 1500, 2000, 2500)){
    print(paste("Processing R =", R, "m =", m, "r =", r))
    result = repeat_experiments(Data, m = m, r = r, rep = 10, log = FALSE)
    df$MSE[k0:(k0+5)] = result$mean
    df$MSE_error[k0:(k0+5)] = result$std
    k0 = k0 + 6
  }
}
write.csv(matrix(df$MSE, 6, 10), file = "time.csv")
