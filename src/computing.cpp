#include "RcppArmadillo.h"
using namespace arma;
#ifdef _OPENMP
#include <omp.h>
#endif

#include <random>
#include <vector>
#include <cstdint>



// [[Rcpp::depends(RcppArmadillo)]]


//' Compute the dense Matrix Multiplication using arma
//'
//' @param A the first matrix
//' @param B the second matrix
//' @return the multiplied matrix
// [[Rcpp::export]]
arma::mat Matmul(const arma::mat& A, const arma::mat& B) {
  arma::mat C = A*B;
  return C;
}

//' Compute the sparse Matrix Multiplication using arma
//'
//' @param A the first matrix, regarded as sparse matrix
//' @param B the second matrix, regarded as dense matrix
//' @return the multiplied matrix
// [[Rcpp::export]]
arma::mat Matmul_spar(const arma::sp_mat& A, const arma::mat& B) {
  arma::mat C = A*B;
  return C;
}


//' Compute the approximation of the kernel matrix
//'
//' @param X 
//' @param Y 
// [[Rcpp::export]]
void Matmul_Approximation(const arma::mat& X, 
                           const arma::mat& Y, 
                           arma::mat& CXX, 
                           arma::mat& CXY, 
                           const arma::mat& Prob, 
                           int r) {
   size_t i, j, k;
   size_t X_rows, X_cols, Y_cols;
   X_rows = X.n_rows;
   X_cols = X.n_cols;
   // Y_rows = Y.n_rows;
   Y_cols = Y.n_cols;
   
   #pragma omp parallel for
   for (i = 0; i < X_cols; ++i) {
     std::random_device rd;
     std::mt19937 gen(rd());
     std::uniform_real_distribution<> dis(0.0, 1.0/r);
     for (j = 0; j < X_rows; ++j){
       if(Prob(j, i) > dis(gen)){
         double xji;
         if(Prob(j, i) > 1.0/r){
           xji = X(j, i);
         }
         else{
           xji = X(j, i) / Prob(j ,i) / r;
         }
         for (k = 0; k < X_cols; ++k) {
           CXX(i, k) += xji * X(j, k);
         }
         for (k = 0; k < Y_cols; ++k) {
           CXY(i, k) += xji * Y(j, k);
         }
       }
     }
   }
 }

// [[Rcpp::export]]
arma::vec RowProb(const arma::mat& X, 
                  const arma::mat& Y, 
                  const std::string& option = "unif") {
  
  arma::vec p;
  
  if (option == "unif") {
    p = arma::ones<arma::vec>(X.n_rows); 
  }
  else if (option == "lev") {
    arma::mat U, V;
    arma::vec s;
    arma::svd_econ(U, s, V, X);
    p = arma::sum(arma::square(U), 1);
  }
  else if (option == "xxopt") {
    p = arma::sum(arma::square(X), 1);
  }
  else if (option == "xyopt") {
    arma::vec a1 = arma::sum(arma::square(X), 1);
    arma::vec a2 = arma::sum(arma::square(Y), 1);
    p = arma::sqrt(a1 % a2);
  }
  else {
    Rcpp::stop("Invalid option. Choose from 'unif', 'lev', 'xxopt', or 'xyopt'.");
  }
  
  p = p / arma::sum(p);
  
  return p;
}


// [[Rcpp::export]]
arma::mat ElementProb(const arma::mat& X, 
                      const arma::mat& Y, 
                      const std::string& option = "xxopt") {
  arma::mat row_norm;
  arma::mat P;
  
  if (option == "xxopt") {
    row_norm = arma::sqrt(arma::sum(arma::square(X), 1));
  } else if (option == "xyopt") {
    row_norm = arma::sqrt(arma::sum(arma::square(Y), 1));
  } else {
    Rcpp::stop("Invalid option. Choose 'xxopt' or 'xyopt'.");
  }

  row_norm = arma::repmat(row_norm, 1, X.n_cols);
  P = arma::abs(X) % row_norm;
  P /= arma::accu(P);
  
  return P;
}




