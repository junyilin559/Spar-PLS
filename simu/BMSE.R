Rcpp::sourceCpp("../src/computing.cpp")
source("../src/SparPLS.R", echo=TRUE)
source("../src/RandomPLS.R", echo=TRUE)
source("../src/IFPLS.R", echo=TRUE)
source("../src/MCDPLS.R", echo=TRUE)

library(Matrix)
library(MASS)
library(pls)
library(RSpectra)

random_orthogonal_matrix = function(dim) {
  Q = qr.Q(qr(matrix(rnorm(dim * dim), nrow = dim)))
  return(Q)
}

Generate_Data = function(n = 50000, p = 100, l = 50, m = 5, R = 0.25){
  P = random_orthogonal_matrix(p)[1:m, ]
  Q = random_orthogonal_matrix(l)[1:m, ]
  beta = t(P) %*% solve(P %*% t(P)) %*% Q
  Sigma = matrix(0, nrow = m, ncol = m)
  for (i in 1:m) {
    for (j in 1:m) {
      Sigma[i, j] = 0.5^(abs(i - j))
    }
  }
  eigen_decomp = eigen(Sigma)
  eigenvalues = eigen_decomp$values
  eigenvectors = eigen_decomp$vectors
  Sigma_1_2 = eigenvectors %*% diag(sqrt(eigenvalues)) %*% t(eigenvectors)
  
  Scale = diag(exp((0:(m-1)) * log(20) / (m-1) / 2))
  SNR_X = R
  SNR_Y = R
  
  T_mat = mvrnorm(n = n, mu = rep(0, m), Sigma = diag(m))
  X0 = T_mat %*% P
  Y0 = T_mat %*% Q
  
  X_std = sd(X0)
  Y_std = sd(Y0)
  
  E_mat = mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p) * (X_std / sqrt(SNR_X))^2)
  F_mat = mvrnorm(n = n, mu = rep(0, l), Sigma = diag(l) * (Y_std / sqrt(SNR_Y))^2)
  
  X = T_mat %*% P + E_mat
  Y = T_mat %*% Q + F_mat
  return(list(X = X, Y = Y, beta = beta, P = P, Q = Q))
}

data_preprocess = function(Data){
  X = Data$X
  p = ncol(X)
  Y = Data$Y
  Prob_spar = ElementProb(X, Y)
  Prob_unif = RowProb(X, Y, option = 'unif')
  Prob_lev = RowProb(X, Y, option = 'lev')
  
  Data$Prob_spar = Prob_spar
  Data$Prob_lev = Prob_lev
  Data$Prob_unif = Prob_unif
  return(Data)
}

get_PLS = function(Data, m = 5, r = 500){
  mse = c(0,0,0,0)
  
  # PLS
  B = PLS.coef(Data$X, Data$Y, m)
  mse[1] = mean((B - Data$beta)*(B-Data$beta))
  
  # RandomPLS by UNIF
  B = RandomPLS(Data$X, Data$Y, m, Data$Prob_unif, r)
  mse[2] = mean((B - Data$beta)*(B-Data$beta))
  
  # RandomPLS by LEV
  B = RandomPLS(Data$X, Data$Y, m, Data$Prob_lev, r)
  mse[3] = mean((B - Data$beta)*(B - Data$beta))
  
  # SparPLS by xxopt
  B = SparPLS(Data$X, Data$Y, m, Data$Prob_spar, r)
  mse[4] = mean((B - Data$beta)*(B - Data$beta))
  
  return(mse)
}

get_PLS_add = function(Data, m = 5, r = 500){
  mse = c(0,0)
  # IFPLS
  B = IFPLS(Data$X, Data$Y, m, r)
  mse[1] = mean((B - Data$beta)*(B-Data$beta))
  # MCDPLS
  B = MCDPLS(Data$X, Data$Y, Data$X, m, r, simu = TRUE)
  mse[2] = mean((B - Data$beta)*(B-Data$beta))
  return(mse)
}

repeat_experiments = function(Data, m = 5, r = 500, 
                              rep = 10, log = TRUE){
  result_mat = matrix(0, rep, 4)
  for(i in c(1:10)){
    result = get_PLS(Data, m, r)
    result_mat[i,] = result
  }
  if(log){
    result_mat = log(result_mat)
  }
  result_mean = apply(result_mat, 2, mean)
  result_mean = c(result_mean[1:3], log(get_PLS_add(Data, m, r)), result_mean[4])
  result_std = apply(result_mat, 2, sd)
  result_std = c(result_std[1:3], 0, 0, result_std[4])
  
  result = list(mean = result_mean, 
                std = result_std)
  return(result)
}

df_init = function(subsample_sizes = c(500, 1000, 1500, 2000, 2500),
                   Method = c("FULL", "UNIF", "LEV", "IF", "MCD", "SPAR"),
                   R_list = c("R1", "R2", "R3"), 
                   M_list = c("M1", "M2", "M3", "M4", "M5")){
  set.seed(123)
  data = expand.grid(
    Method = Method,
    subsample = subsample_sizes,
    R = R_list,
    M = M_list
  )
  
  data$MSE = rnorm(nrow(data), mean = 10, sd = 2)
  data$MSE_error = rnorm(nrow(data), mean = 1, sd = 0.5)
  data$subplot_label = paste(data$R, data$M, sep = ",")
  return(data)
}

set.seed(1234)
df = df_init()
k0 = 1
for(R in c(0.25, 0.5, 1)){
  for(m in c(5, 10, 15, 20, 25)){
    Data = Generate_Data(m = m, R = R)
    Data = data_preprocess(Data)
    for(r in c(500, 1000, 1500, 2000, 2500)){
      print(paste("Processing R =", R, "m =", m, "r =", r))
      result = repeat_experiments(Data, m = m, r = r, rep = 10, log = TRUE)
      df$MSE[k0:(k0+5)] = result$mean
      df$MSE_error[k0:(k0+5)] = result$std
      k0 = k0 + 6
    }
  }
}



# Plot
library(ggplot2)
p <- ggplot(df, aes(x = subsample, y = MSE, group = Method, linetype = Method, color = Method)) +
  geom_line(size = 0.8) +  
  geom_point(aes(shape = Method), size = 2) +
  geom_errorbar(aes(ymin = MSE - MSE_error, ymax = MSE + MSE_error), 
                width = 80, size = 0.8, linetype = "solid", show.legend = FALSE) +  
  facet_wrap(~subplot_label, nrow = 3, ncol = 5, scales = "free_y", labeller = label_value) +  
  labs(x = "Number of Subsamples", y = "log(MSE)") +  
  scale_color_manual(
    values = c("FULL" = "black", "SPAR" = "red", "UNIF" = "grey", "LEV" = "blue", "IF" = "#ae017e", MCD = "#a6cee3"), 
    name = "Method"  
  ) +
  theme_minimal() +  
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(1, 1, 1, 1),
    #legend.title = element_blank(),
    legend.title = element_text(size = 14),  
    #legend.background = element_rect(fill = alpha("white", 0.6), color = "black"), 
    legend.key.width = unit(2, "line"), 
    legend.key.height = unit(1, "line"), 
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 10),  
    axis.ticks = element_line(color = "black", size = 0.5),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    plot.margin = margin(10, 10, 10, 10)  
  ) + 
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         linetype = guide_legend(nrow = 1, byrow = TRUE),
         shape = guide_legend(nrow = 1, byrow = TRUE))
p


