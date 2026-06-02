IF_tor = function(X, Y, T_, B){
  CXX = t(X) %*% X
  lev = colSums(solve(CXX, t(X)) * t(X))
  Res = Y - X %*% B
  CTY = t(T_) %*% Res
  R_square = colSums(T_ %*% CTY * Res)
  tor = 2 / (1 - lev) * R_square
  return(tor)
}


IFPLS = function(X, Y, m, n_sample){
  X0 = scale(X)
  Y0 = scale(Y)
  model <- plsr(Y ~ X, ncomp = m, validation = "none", center = F)
  B = coef(model, ncomp = m, intercept = FALSE)[,,1]
  T_ = model$scores
  IFTOR = IF_tor(X, Y, T_, B)
  index = order(IFTOR, decreasing = TRUE)[1:n_sample]
  X_sub = X[index,]
  Y_sub = Y[index,]
  
  return(PLS.coef(X_sub, Y_sub, m))
}