#include <math.h>
#include <RcppArmadillo.h>
#include <mvt.h>
#include <mvnorm.h>
#include <ll.h> // straight copy from gmvjoint

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;

// Testing RcppDist stuff
// [[Rcpp::export]]
void test(uword n, vec& x, vec& mu, mat& S){
  Rcout << "n, x, mu, S:" << n << std::endl;
  x.print();
  mu.print();
  S.print();
  
  Rcout << "MVN from RcppDist" << std::endl;
  // Transform x to a matrix 
  mat xx = reshape(x, 1, x.size());
  xx.print();
  vec dmv = dmvnorm(xx, mu, S);
  dmv.print();
  double dmv2 = as_scalar(dmvnorm(xx, mu, S));
  Rcout << "dmv2 = " << dmv2 << std::endl;
  
  // Random draw -->
  Rcout << "Random draw ->" << std::endl;
  mat rmvn = rmvnorm(n, mu, S);
  rmvn.print();
  vec vrmvn = vectorise(rmvn, 0);
  vrmvn.print();
  
  Rcout << "----" << std::endl;
  Rcout << "MVt from RcppDist" << std::endl;
  double ddmvt = as_scalar(dmvt(xx, mu, S, 4.));
  Rcout << "dmvt = " << ddmvt << std::endl;
  
  // Random draw ->
  vec vrmvt = vectorise(rmvt(n, mu, S, 4.), 0);
  Rcout << "Random draw ->" << std::endl;
  vrmvt.print();
}

// Importing from gmvjoint ------------------------------------------------
// I'm too thick to work out how to declare/"borrow" directly!
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

double logfti(const arma::vec& b, const arma::rowvec& S, const arma::mat& SS, const arma::rowvec& Fi, const arma::mat& Fu,
              const double l0i, const arma::rowvec& haz, const int Delta, const arma::vec& gamma_rep, const arma::vec& zeta){
  double temp = 0.0;
  if(Delta == 1) temp = log(l0i);
  return as_scalar(
    temp + Delta * (S * zeta + Fi * (gamma_rep % b)) - haz * exp(SS * zeta + Fu * (gamma_rep % b))
  );
}
  
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

// Functions to reshape vec -> mat, and simulate from either normal or t distn.
// [[Rcpp::export]]
vec draw_b(const vec& b_current, const mat& Sigmatuned, const std::string& b_dist, const double df){
  uword n = 1;
  if(b_dist == "normal"){ // Normal, or...
    return vectorise(rmvnorm(n, b_current, Sigmatuned));
  }else{                  // multivariate t draw
    return vectorise(rmvt(n, b_current, Sigmatuned, df), 0);
  }
}

double eval_b(vec& b, const vec& mu, const mat& Sigmatuned, const std::string& b_dist, const double df){
  // Transform to a matrix, since RcppDist operates under this assumption...
  mat bb = reshape(b, 1, b.size());
  if(b_dist == "normal"){ // Normal, or...
    return as_scalar(dmvnorm(bb, mu, Sigmatuned, true));
  }else{                  // multivariate t draw
    return as_scalar(dmvt(bb, mu, Sigmatuned, df, true));
  }
}

// Metropolis algorithm for sampling from f(b_i|...; Omega) ---------------
// [[Rcpp::export]]
List Metropolis_Hastings(const arma::vec& b, const List& Omega,
                         const List& Y, const List& X, const List& Z, const List& W,
                         const List& family, const int Delta, const arma::rowvec& S, const arma::rowvec& Fi,
                         const double l0i, const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz,
                         const arma::vec& gamma_rep, const List& beta_inds, const List& b_inds, const arma::uword K,
                         const arma::uword q, const int burnin, const int N, const arma::mat& Sigma,
                         const std::string& b_dist, const double df,
                         const double tune){
  // Unpack Omega
  mat D = Omega["D"];
  List sigma = Omega["sigma"];
  vec beta = Omega["beta"];
  vec zeta = Omega["zeta"];
  mat Sigmatuned = Sigma * tune;
  // MC stuff
  int iters = burnin + N;
  mat out = mat(q, iters);
  out.col(0) = b;
  // Start
  int j = 1, num_accepts = 0;
  while(j < iters){
    Rcpp::checkUserInterrupt(); // So function is interruptable; might be useful.
    double U = randu();
    vec b_current = out.col(j - 1);
    // Draw from proposal distribution (coding _both_ )
    vec b_proposal = draw_b(b, Sigmatuned, b_dist, df);
    // MH part, log(f|...)
    double logf_current = -1. * joint_density(b_current, Y, X, Z, W, beta, D, sigma, family, Delta,
                                              S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double logf_proposal = -1. * joint_density(b_proposal, Y, X, Z, W, beta, D, sigma, family, Delta,
                                              S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double log_firstbit = logf_proposal - logf_current;
    // Underlying density
    double dens_current = eval_b(b_current,  b, Sigmatuned, b_dist, df);
    double dens_proposal = eval_b(b_proposal, b, Sigmatuned, b_dist, df);
    double log_secondbit = dens_current - dens_proposal;
    double P = std::min(exp(log_firstbit - log_secondbit), 1.);
    if(U < P){
      b_current = b_proposal;
      if(j > burnin) num_accepts ++;
    }
    out.col(j) = b_current;
    j++;
  }

  // Remove burnin -----------
  out.shed_cols(0, burnin - 1);
  return List::create(_["walks"] = out,
                      _["burnin"] = burnin,
                      _["N"] = N,
                      _["AcceptanceRate"] = (double)num_accepts/(double)N);
}
