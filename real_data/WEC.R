# Download the Wave Energy Converters dataset from:  
# https://archive.ics.uci.edu/dataset/494/wave+energy+converters  
# Ensure the dataset is stored in the './data' directory.  

Rcpp::sourceCpp("../src/computing.cpp")
source("../src/SparPLS.R", echo=TRUE)
source("../src/RandomPLS.R", echo=TRUE)
source("../src/IFPLS.R", echo=TRUE)
source("../src/MCDPLS.R", echo=TRUE)

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
  mse = c(0, 0, 0, 0, 0, 0)
  
  X_train = Data$X_train
  Y_train = Data$Y_train
  X_test = Data$X_test
  Y_test = Data$Y_test
  
  # PLS
  B <- PLS.coef(X_train, Y_train, m)
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[1] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[1] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  
  # RandomPLS by UNIF
  B <- RandomPLS(X_train, Y_train, m, Data$Prob_unif, r)
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[2] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[2] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  
  # RandomPLS by LEV
  B <- RandomPLS(X_train, Y_train, m, Data$Prob_lev, r)
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[3] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[3] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  
  # IFPLS
  B = IFPLS(X_train, Y_train, m, r)
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[4] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[4] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  
  # MCDPLS
  #B = MCDPLS(X_train, Y_train, X_test, m, r)
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[5] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[5] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  
  # SparPLS by xxopt
  dense_rate = mean(X_train != 0)
  B = SparPLS(X_train, Y_train, m, Data$Prob_spar, round(r*dense_rate))
  Y_pred = X_test %*% B
  if(task == 'Regression'){
    mse[6] = mean((Y_pred - Y_test)^2) 
  }else{
    mse[6] = mean(apply(Y_pred, 1, which.max) == apply(Y_test, 1, which.max)) 
  }
  return(mse)
}

cross_validate = function(Data, m = 5, r = 500, 
                          cv_folds = 10, log = TRUE, 
                          task = 'Regression'){
  mse_fold = matrix(0, cv_folds, 6)
  
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
    
    mse_fold[fold, ] = result
  }
  
  if(log){
    mse_fold = log(mse_fold)
  }
  
  result_mean = apply(mse_fold, 2, mean)
  result_std = apply(mse_fold, 2, sd)
  
  result = list(mean = result_mean, std = result_std)
  return(result)
}


set.seed(1234)
results_list <- list()
for(loc in c('Perth', 'Sydney')){
  for(num in c(49, 100)){
    Data = get_Data(loc = loc, num = num)
    n = nrow(Data$X)
    r0 = round(n / 50)
    for(m in c(5, 10, 15, 20)){
      for(r in c(r0, 2*r0, 3*r0, 4*r0, 5*r0)){
        print(paste('loc: ', loc, '; num: ', num, '; m: ', m, '; r: ', r))
        result = cross_validate(Data, m = m, r = r, 
                                cv_folds = 10, log = TRUE, 
                                task = 'Regression')
        print(result)
        
        methods <- c("FULL", "UNIF", "LEV", "IF", "MCD", "SPAR")
        for (i in seq_along(methods)) {
          results_list <- append(results_list, list(data.frame(
            loc = loc,
            num = num,
            m = m,
            r = r,
            method = methods[i],
            mean = result$mean[i],
            std = result$std[i]
          )))
        }
      }
    }
  }
}

final_df <- do.call(rbind, results_list)

saveRDS(final_df, "WEC1.rds")


final_df = readRDS("WEC1.rds")
final_df$r = rep(c(1:5), each = 6, times = 16)
final_df$subplot_label <- paste0("Loc: ", final_df$loc, ", Num: ", final_df$num)
final_df$subplot_label = rep(c('C1', 'C2','C3','C4'), each = 120)
final_df$m = rep(c('M1', 'M2','M3','M4'), each = 30,times = 4)
final_df$sub_title = paste0(final_df$subplot_label, ", ", final_df$m)

final_df$method <- factor(final_df$method,
                          levels = c("FULL", "UNIF", "LEV", "IF", "MCD", "SPAR")
)
final_df <- final_df |>
  dplyr::filter(!(method == "MCD" & subplot_label %in% c("C1", "C3")))


# Plot
library(ggplot2)
p <- ggplot(final_df, aes(x = r, y = mean, group = method, linetype = method, color = method)) +
  geom_line(size = 1) +  
  geom_point(aes(shape = method), size = 2) +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), 
                width = 0.2, size = 1, linetype = "solid", show.legend = FALSE) +  
  facet_wrap(~sub_title, nrow = 4, ncol = 4, scales = "free_y", shrink = TRUE, labeller = label_value) +  
  labs(x = "Number of Subsamples", 
       y = "log(PMSE)",  
       color = "method") +  
  scale_color_manual(
    values = c("FULL" = "black", "SPAR" = "red", "UNIF" = "grey", "LEV" = "blue", "IF" = "#ae017e", MCD = "#a6cee3"), 
    name = "method"  
  ) +
  scale_x_continuous(
    breaks = function(x) seq(1,5),  
    labels = function(x) paste0(c('', 2, 3, 4, 5), "r")  
  ) +
  theme_minimal() + 
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(1, 1, 1, 1),
    legend.box.margin = margin(1, 1, 1, 1),
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






