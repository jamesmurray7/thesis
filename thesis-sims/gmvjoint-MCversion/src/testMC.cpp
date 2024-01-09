#include <math.h>
#include <RcppArmadillo.h>
#include "LOGLIK.h"
#include "COMPDATASCORE.h"
#include "BETAUPDATES.h"

// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using namespace arma;

// Create etas for subject i
// [[Rcpp::export]]
List make_eta(const List& X, const List& Z, const arma::vec& beta, const arma::vec& b,
              const List& beta_inds, const List& b_inds){
  uword K = X.size();
  List out(K);
  for(uword k = 0; k < K; k++){
    uvec beta_inds_k = beta_inds[k], b_inds_k = b_inds[k];
    mat X_k = X[k], Z_k = Z[k];
    vec beta_k = beta.elem(beta_inds_k), b_k = b.elem(b_inds_k);
    out[k] = X_k * beta_k + Z_k * b_k;
  }
  return out;
}

// [[Rcpp::export]]
double logfti(const arma::vec& b, const arma::rowvec& S, const arma::mat& SS, const arma::rowvec& Fi, const arma::mat& Fu,
             const double l0i, const arma::rowvec& haz, const int Delta, const arma::vec& gamma_rep, const arma::vec& zeta){
 double temp = 0.0;
 if(Delta == 1) temp = log(l0i);
 
 return as_scalar(
   temp + Delta * (S * zeta + Fi * (gamma_rep % b)) - haz * exp(SS * zeta + Fu * (gamma_rep % b))
   );
}

// [[Rcpp::export]]
double joint_density(const arma::vec& b, 
                     const List& Y, const List& X, const List& Z, const List& W,
                     const arma::vec& beta, const arma::mat& D, const List& sigma, const List& family,
                     const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, const double l0i,
                     const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                     const arma::vec& gamma_rep, const arma::vec& zeta,
                     const List& beta_inds, const List& b_inds, const arma::uword K){
  double ll = 0.0; // Compile longitudinal parts ---------
  List eta = make_eta(X, Z, beta, b, beta_inds, b_inds);
  for(uword k = 0; k < K; k++){
    vec Yk = Y[k];
    std::string f = family[k];
    vec sigmak = sigma[k];
    vec eta_k = eta[k];
    if(f == "gaussian"){
      ll += ll_Gaussian(eta_k, Yk, as_scalar(sigmak));
    }else if(f == "binomial"){
      ll += ll_Binomial(trunc_exp(eta_k)/(1. + trunc_exp(eta_k)), Yk);
    }else if(f == "poisson"){
      ll += ll_Poisson(trunc_exp(eta_k), Yk);
    }else if(f == "negbin"){
      mat W_k = W[k];
      ll += ll_NegBin(trunc_exp(eta_k), Yk, trunc_exp(W_k * sigmak));
    }else if(f == "genpois"){
      mat W_k = W[k];
      ll += ll_GenPois(trunc_exp(eta_k), Yk, W_k * sigmak);
    }else if(f == "Gamma"){
      mat W_k = W[k];
      ll += ll_Gamma(trunc_exp(eta_k), Yk, trunc_exp(W_k * sigmak)); 
    }
  }
  // uword q = b.size();
  int q = b.size();
  rowvec zz = zeros<rowvec>(q);
  
  double ll_Ti = logfti(b, S, SS, Fi, Fu, l0i, haz, Delta, gamma_rep, zeta);
  // double ll_RE = as_scalar(dmvnrm_arma_fast(b.t(), zz, D, true));
  double ll_RE = as_scalar(-(double)q/2. * log2pi - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b);
  
  return -1.0 * (ll + ll_RE + ll_Ti);
}

// [[Rcpp::export]]
arma::vec joint_density_ddb(const arma::vec& b, 
                            const List& Y, const List& X, const List& Z, const List& W,
                            const arma::vec& beta, const arma::mat& D, const List& sigma, const List& family,
                            const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, const double l0i,
                            const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                            const arma::vec& gamma_rep, const arma::vec& zeta,
                            const List& beta_inds, const List& b_inds, const arma::uword K){
  uword q = b.size();
  vec Score(q);
  List eta = make_eta(X, Z, beta, b, beta_inds, b_inds);
  for(uword k = 0; k < K; k++){
    vec Yk = Y[k];
    std::string f = family[k];
    vec sigmak = sigma[k];
    vec eta_k = eta[k];
    mat Z_k = Z[k];
    uvec b_inds_k = b_inds[k];
    if(f == "gaussian"){
      Score.elem(b_inds_k) += Z_k.t() * score_Gaussian(eta_k, Yk, as_scalar(sigmak));
    }else if(f == "binomial"){
      Score.elem(b_inds_k) += Z_k.t() * score_Binomial(eta_k, Yk);
    }else if(f == "poisson"){
      Score.elem(b_inds_k) += Z_k.t() * score_Poisson(eta_k, Yk);
    }else if(f == "negbin"){
      mat W_k = W[k];
      Score.elem(b_inds_k) += Z_k.t() * score_NegBin(eta_k, Yk, trunc_exp(W_k * sigmak));
    }else if(f == "genpois"){
      mat W_k = W[k];
      Score.elem(b_inds_k) += Z_k.t() * score_GenPois(eta_k, Yk, W_k * sigmak);
    }else if(f == "Gamma"){
      mat W_k = W[k];
      Score.elem(b_inds_k) += Z_k.t() * score_Gamma(eta_k, Yk, trunc_exp(W_k * sigmak));
    }
  }
  
  vec Score_Ti_b = Score_Ti(b, S, SS, Fi, Fu, l0i, haz, Delta, gamma_rep, zeta);
  vec Score_b = grad_b(b, D);
  
  return -1.0 * (Score + Score_b + Score_Ti_b);
}
// This isn't actually needed!
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
List make_binExp(const vec& Y, const mat& X, const mat& Z, 
            const vec& beta, const mat& draws){
  uword N = draws.n_rows;
  vec out1(Y.size()), out2(Y.size());
  for(uword i = 0; i < N; i++){
    vec b_i = draws.row(i).t();
    vec yyy = exp(X * beta + Z * b_i);
    out1 += yyy / (1. + yyy);
    out2 += yyy / (pow(1. + yyy, 2.));
  }
  return List::create(_["dot"] = out1 / N,
                      _["ddot"] = out2/ N);
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

// [[Rcpp::export]]
double Egammazeta(const arma::vec& gammazeta, const arma::mat& draws, const vec& Eb,
                  const arma::rowvec& S, const arma::mat& SS, const arma::mat& Fu, 
                  const arma::rowvec& Fi, const arma::vec& haz, const int Delta, 
                  const List& inds, const arma::uword K){
  uword q = Fu.n_cols;
  vec g = gammazeta.head(K), z = gammazeta.subvec(K, gammazeta.size() - 1), gammas(q);
  for(uword k = 0; k < K; k++){
    double gk = g[k];
    uvec inds_k = inds[k];
    gammas.elem(inds_k) += gk;
  }
  vec Esurvpart = make_Esurvpart(draws, SS, Fu, gammas, z);
  double rhs = as_scalar(haz.t() * Esurvpart);
  return as_scalar(Delta * (S * z + Fi * (Eb % gammas))) - rhs;
}

// [[Rcpp::export]]
arma::vec Sgammazeta(arma::vec& gammazeta, const arma::mat& draws, const vec& Eb,
                     arma::rowvec& S, arma::mat& SS, arma::mat& Fu, arma::rowvec& Fi,
                     arma::vec& haz, int Delta, Rcpp::List b_inds, uword K, double eps){
  int ps = gammazeta.size();
  vec out = vec(ps), epsvec = vec(ps);
  for(int i = 0; i < ps; i++){
    epsvec[i] = eps;
    double xi = 1.;
    vec ge1 = gammazeta + xi * epsvec, ge2 = gammazeta - xi * epsvec;
    double f1 = Egammazeta(ge1, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K),
           f2 = Egammazeta(ge2, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K);
    out[i] = (f1 - f2) / (2 * eps);
    epsvec[i] = 0.;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat Hgammazeta(arma::vec& gammazeta, const arma::mat& draws, const vec& Eb,
                     arma::rowvec& S, arma::mat& SS, arma::mat& Fu, arma::rowvec& Fi,
                     arma::vec& haz, int Delta, Rcpp::List b_inds, uword K, double eps){
  int ps = gammazeta.size();
  mat out = zeros<mat>(ps, ps);
  vec epsvec = vec(ps, fill::value(eps));
  mat epsmat = diagmat(epsvec);
  double eps2 = pow(eps, 2.), eps42 = 4. * pow(eps, 2.);
  double f0 = Egammazeta(gammazeta, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K);
  for(int i = 0; i < (ps - 1); i++){
    vec eps_i = epsmat.col(i);
    // Diagonal terms
    vec ge_i1 = gammazeta + eps_i, ge_i2 = gammazeta - eps_i;
    double f1 = Egammazeta(ge_i1, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K),
           f2 = Egammazeta(ge_i2, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K);
    out(i,i) = (f1 - 2. * f0 + f2)/eps2;
    // Off-diagonal
    for(int j = (i + 1); j < ps; j++){
      vec eps_j = epsmat.col(j);
      vec ge_ij1 = gammazeta + eps_i + eps_j,
          ge_ij2 = gammazeta + eps_i - eps_j,
          ge_ij3 = gammazeta - eps_i + eps_j,
          ge_ij4 = gammazeta - eps_i - eps_j;
      double f1 = Egammazeta(ge_ij1, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K),
             f2 = Egammazeta(ge_ij2, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K),
             f3 = Egammazeta(ge_ij3, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K),
             f4 = Egammazeta(ge_ij4, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K);
      out(i,j) = (f1 - f2 - f3 + f4)/eps42;
      out(j,i) = out(i,j);
    }
  }
  // Calculate ps, psth item
  int last = ps-1;
  vec eps_i = epsmat.col(ps-1);
  vec ge_i1 = gammazeta + eps_i, ge_i2 = gammazeta - eps_i;
  double f1 = Egammazeta(ge_i1, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K),
         f2 = Egammazeta(ge_i2, draws, Eb, S, SS, Fu, Fi, haz, Delta, b_inds, K);
  out(last, last) = (f1 - 2. * f0 + f2)/eps2;
  return out;
}

// [[Rcpp::export]]
arma::vec lambda_update(Rcpp::List draws, Rcpp::List Fu, Rcpp::List SS,
                        Rcpp::List survtimes, arma::vec& gamma_rep, arma::vec& zeta, arma::vec& nev){
  uword n = survtimes.size(), unique_times = nev.size();
  vec out = zeros<arma::vec>(unique_times);
  
  for(uword i = 0; i < n; i++){
    
    vec st_i = survtimes[i];
    mat SS_i = SS[i], Fu_i = Fu[i], draws_i = draws[i];
    
    int num_survived = st_i.size();
    if(num_survived == 0) continue;
    
    vec Esurvpart = make_Esurvpart(draws_i, SS_i, Fu_i, gamma_rep, zeta);
    
    int ui = Fu_i.n_rows;
    out.subvec(0, ui - 1) += Esurvpart;
    
  }
  return nev/out;
}