PLS.coef = function(X, Y, m){
  model <- plsr(Y ~ X, ncomp = m, validation = "none", center = F)
  B = coef(model, ncomp = m, intercept = FALSE)[,,1]
  return(B)
}

RandomPLS = function(X, Y, m, Prob, r){
  #Row-wise PLS by kernel method
  
  # Initialize
  n = nrow(X)
  p = ncol(X)
  l = ncol(Y)
  sampled_indices <- sample(1:n, prob = Prob, size = r, replace = TRUE)
  X_sub = X[sampled_indices, , drop = FALSE] / sqrt(Prob[sampled_indices])
  Y_sub = Y[sampled_indices, , drop = FALSE] / sqrt(Prob[sampled_indices])
  B = PLS.coef(X_sub, Y_sub, m)
  return(B)
}
