C_alpha <- function(p, alpha) {
  q <- qchisq(alpha, df = p)
  alpha / pchisq(q, df = p + 2)
}

MCDPLS = function(X_train, Y_train, X_test, m, r, simu = FALSE){
  X = rbind(X_train, X_test)
  n = nrow(X)
  n_train = nrow(X_train)
  n_test = nrow(X_test)
  p = ncol(X)
  if(simu){
    n = n_test
    n_train = 0
    X = X_test
  }
  
  
  # loop
  converged = FALSE
  h = p+10
  while (!converged) {
    C_step_convergence = FALSE
    H = X_test
    while(!C_step_convergence){
      T0 = colMeans(H)
      C_al = C_alpha(p, h/n)
      S =  C_al * cov(H)
      
      MD = sqrt(colSums((ginv(S) %*% t(X - T0)) * t(X - T0)))
      H_new = X[order(MD)[1:h], ]
      
      if (abs(det(C_al * cov(H_new)) - det(S)) < 1e-12) {
        C_step_convergence = TRUE
      } else {
        H = H_new
      }
    }
    
    T0 = colMeans(H)
    C_al = C_alpha(p, h/n)
    S =  C_al * cov(H)
    MD = sqrt(colSums((ginv(S) %*% t(X - T0)) * t(X - T0)))
    if(all(MD[(n_train+1):n] <= qchisq(0.975, p)) | (h+10 > n)){
      converged = TRUE
    } else {
      h = h + 10
    }
  }
  
  
  #MD_train = sqrt(colSums(solve(S2, t(X_train - T0)) * t(X_train - T0)))
  MD_train = sqrt(colSums((ginv(S) %*% t(X_train - T0)) * t(X_train - T0)))
  #X_sub = X_train[which(MD_train < qchisq(0.975, p)),]
  #Y_sub = Y_train[which(MD_train < qchisq(0.975, p)),]
  X_sub = X_train[order(MD_train)[1:r],]
  Y_sub = Y_train[order(MD_train)[1:r],]
  return(PLS.coef(X_sub, Y_sub, m))
}
