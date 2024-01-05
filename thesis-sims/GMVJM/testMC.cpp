#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat sum_tcrossprods(const mat& X) {
  uword N = X.n_rows, q = X.n_cols;
  
  mat out(q, q, fill::zeros);
  
  // Iterate over each row of the matrix
  for (uword i = 0; i < N; i++) {
    // Extract the i-th row as a column vector
    vec x_i = X.row(i).t();
    
    // Compute the outer product of the row vector with itself
    mat xxT = x_i * x_i.t();
    
    // Add the outer product to the result matrix
    out += xxT;
  }
  
  return out;
}

// [[Rcpp::export]]
List make_Expeta(const List& X, const List& Z, const arma::vec& beta, 
                 const arma::mat& draws,
                 const List& beta_inds, const List& b_inds){
  uword K = X.size(), N = draws.n_rows, q = draws.n_cols;
  List out(K);
  for(uword k = 0; k < K; k++){
    uvec beta_inds_k = beta_inds[k], b_inds_k = b_inds[k];
    mat X_k = X[k], Z_k = Z[k];
    vec beta_k = beta.elem(beta_inds_k);
    vec this_out = vec(X_k.n_rows);
    for(uword i = 0; i < N; i++){
      vec b_i = draws.row(i).t();
      vec b_k = b_i.elem(b_inds_k);
      this_out += exp(X_k * beta_k + Z_k * b_k);
    }
    out[k] = this_out / N;
  }
  return out;
}

// [[Rcpp::export]]
double Eymeta(const vec& Y, const mat& X, const mat& Z, 
              const vec& beta, const mat& draws){
  uword N = draws.n_rows, q = draws.n_cols;
  double out = 0.;
  for(uword i = 0; i < N ; i++){
    vec b_i = draws.row(i).t();
    vec yyy = Y - (X * beta + Z * b_i);
    out += as_scalar(yyy.t() * yyy);
  }
  return out / N;
}

// [[Rcpp::export]]
mat EymetaT(const vec& Y, const mat& X, const mat& Z, 
            const vec& beta, const mat& draws){
  uword N = draws.n_rows, P = beta.size();
  mat out(P, P);
  for(uword i = 0; i < N ; i++){
    vec b_i = draws.row(i).t();
    vec yyy = Y - (X * beta + Z * b_i);
    out += yyy * yyy.t();
  }
  return out / N;
}

// [[Rcpp::export]]
vec make_Esurvpart(const mat& draws,
                   const mat& SS, const mat& Fu, 
                   const vec& gamma_rep, const vec& zeta){
  uword N = draws.n_rows, u = Fu.n_rows, q = Fu.n_cols;
  vec out(u), Sz = SS * zeta;
  mat gmat = diagmat(gamma_rep);
  mat Q = Fu * gmat;
  for(int i = 0; i < N; i++){
    vec b_i = draws.row(i).t();
    out += exp(Sz + Q * b_i);
  }
  return out/N;
}