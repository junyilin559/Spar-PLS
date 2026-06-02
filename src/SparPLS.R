SparPLS = function(X, Y, m, Prob, num_sample){
  # Spar-PLS by kernel method
  
  # Initialize
  n = nrow(X)
  p = ncol(X)
  l = ncol(Y)
  
  CXX = matrix(0, p, p)
  CXY = matrix(0, p, l)
  Matmul_Approximation(X, Y, CXX, CXY, Prob, num_sample*p)
  CYX = t(CXY)
  chg = diag(p)
  R = matrix(0, p, m)
  Q = matrix(0, l, m)
  
  # loop
  for(i in c(1:m)){
    svd_result <- svds(CYX, k = 1)
    #svd_result <- svd(CYX)
    v = svd_result$v[,1]
    r = chg %*% v
    R[, i] = r
    normt = (t(r) %*% CXX %*% r)
    p_vec = CXX %*% r / normt[1]
    p_vec1 = t(chg) %*% t(CXX) %*% r / normt[1]
    q_vec = CYX %*% r / normt[1]
    
    Q[, i] = q_vec
    chg = chg - chg %*% matrix(v, ncol=1) %*% t(p_vec)
    #chg = chg %*% (diag(p) - matrix(v) %*% t(p_vec)) 
    CYX = CYX - q_vec %*% t(p_vec) * normt[1]
  }
  
  B = R %*% t(Q)
  return(B)
}
