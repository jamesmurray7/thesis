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

// Create taus for subject i
// [[Rcpp::export]]
List make_tau(const List& Z, const List& Sigma){
  uword K = Sigma.size();
  List out(K);
  for(uword k = 0; k < K; k++){
    mat Zk = Z[k], Sk = Sigma[k];
    out[k] = sqrt(diagvec(Zk * Sk * Zk.t()));
  }
  return out;
}


//' @keywords internal
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

// Update to beta ---------------------------------------------------------
// [[Rcpp::export]]
arma::vec Sbeta(const arma::vec& beta, const List& X, const List& Y, const List& Z, const List& W,
                const arma::vec& b, const List& sigma, const List& family, 
                const List& beta_inds, const List& b_inds, const arma::uword& K,
                const List& tau, const arma::vec& w, const arma::vec& v){
  uword P = beta.size();
  vec Score(P);
  List eta = make_eta(X, Z, beta, b, beta_inds, b_inds);
  for(uword k = 0; k < K; k++){
    vec Yk = Y[k], sigmak = sigma[k], tauk = tau[k], etak = eta[k];
    mat Wk = W[k], Xk = X[k];
    uvec k_inds = beta_inds[k];
    std::string ff = family[k];
    if(ff == "gaussian"){
      Score.elem(k_inds) += Xk.t() * d_eta_Gaussian(etak, Yk, as_scalar(sigmak));
    }else if(ff == "poisson"){
      Score.elem(k_inds) += Xk.t() * d_eta_Poisson(etak, Yk);
    }else if(ff == "binomial"){
      Score.elem(k_inds) += Xk.t() * d_eta_Binomial(etak, Yk, tauk, w, v);
    }else if(ff == "negbin"){
      vec phi = trunc_exp(Wk * sigmak);
      Score.elem(k_inds) += Xk.t() * d_eta_NegBin(etak, Yk, phi, tauk, w, v);
    }else if(ff == "Gamma"){
      vec shape = trunc_exp(Wk * sigmak);
      Score.elem(k_inds) += Xk.t() * d_eta_Gamma(etak, Yk, shape);
    }else if(ff == "genpois"){
      vec phi = Wk * sigmak;
      Score.elem(k_inds) += Xk.t() * d_eta_GenPois(etak, Yk, phi, tauk, w, v);
    }else{
      Rcout << "Unknown family " << ff << std::endl;
      continue;
    }
  }
  return Score;
}

// [[Rcpp::export]]
arma::mat Hbeta(const arma::vec& beta, const List& X, const List& Y, const List& Z, const List& W,
                const arma::vec& b, const List& sigma, const List& family, 
                const List& beta_inds, const List& b_inds, const arma::uword& K,
                const List& tau, const arma::vec& w, const arma::vec& v){
  uword P = beta.size();
  mat Hessian(P, P);
  List eta = make_eta(X, Z, beta, b, beta_inds, b_inds);
  for(uword k = 0; k < K; k++){
    vec Yk = Y[k], sigmak = sigma[k], tauk = tau[k], etak = eta[k];
    mat Wk = W[k], Xk = X[k];
    uvec k_inds = beta_inds[k];
    std::string ff = family[k];
    if(ff == "gaussian"){
      Hessian.submat(k_inds, k_inds) = Xk.t() * d2_eta_Gaussian(Yk, as_scalar(sigmak)) * Xk;
    }else if(ff == "poisson"){
      Hessian.submat(k_inds, k_inds) = form_hess(d2_eta_Poisson(etak, Yk), Xk);
    }else if(ff == "binomial"){
      Hessian.submat(k_inds, k_inds) = form_hess(d2_eta_Binomial(etak, tauk, w, v), Xk);
    }else if(ff == "negbin"){
      vec phi = trunc_exp(Wk * sigmak);
      Hessian.submat(k_inds, k_inds) = form_hess(d2_eta_NegBin(etak, Yk, phi, tauk, w, v), Xk);
    }else if(ff == "Gamma"){
      vec shape = trunc_exp(Wk * sigmak);
      Hessian.submat(k_inds, k_inds) = form_hess(d2_eta_Gamma(etak, Yk, shape), Xk);
    }else if(ff == "genpois"){
      vec phi = Wk * sigmak;
      Hessian.submat(k_inds, k_inds) = form_hess(d2_eta_GenPois(etak, Yk, phi, tauk, w, v), Xk);
    }else{
      Rcout << "Unknown family " << ff << std::endl;
      continue;
    }
  }
  return Hessian;
}

// Update to shape/dispersion parameter sigma -----------------------------

// Take three by central differencing (via pracma::grad and pracma::hessian)
// Gamma, shape = exp{W * sigma}
// [[Rcpp::export]]
double appxE_Gammasigma(const arma::vec& sigma, const arma::vec& eta, const arma::vec& Y, 
                        const arma::vec& tau, const arma::mat& W,
                        const arma::vec& w, const arma::vec& v){
  uword gh = w.size();
  vec shape = trunc_exp(W * sigma);
  double out = 0.;
  for(uword l = 0; l < gh; l++){
    vec this_mu = trunc_exp(eta + tau * v[l]);
    out += w[l] * ll_Gamma(this_mu, Y, shape);
  }
  return out;
}

// NegBin, (over)dispersion = exp{W * sigma}
// [[Rcpp::export]]
double appxE_NegBinsigma(const arma::vec& sigma, const arma::vec& eta, const arma::vec& Y, 
                         const arma::vec& tau, const arma::mat& W,
                         const arma::vec& w, const arma::vec& v){
  uword gh = w.size(), m = eta.size();
  vec phi = trunc_exp(W * sigma), Exp(m);
  vec p1 = lgamma(Y + phi) - lgamma(phi) - lgamma(Y + 1.) + phi % log(phi);
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    Exp += w[l] * log(trunc_exp(this_eta) + phi);
  }
  return sum(p1 - (phi + Y) % Exp);
}

// GenPois, dispersion = W * sigma
// [[Rcpp::export]]
double appxE_GenPoissigma(const arma::vec& sigma, const arma::vec& eta, const arma::vec& Y, 
                          const arma::vec& tau, const arma::mat& W,
                          const arma::vec& w, const arma::vec& v){
  // uword gh = w.size(), m = eta.size();
  // vec phi = W * sigma, Exp(m);
  // vec lfac = lgamma(Y + 1.), 
  //     frac = (trunc_exp(eta) + phi % Y) / (1. + phi);
  // for(uword l = 0; l < gh; l++){
  //   vec this_eta = eta + tau * v[l];
  //   Exp += w[l] * log(trunc_exp(this_eta) + phi % Y);
  // }
  // return sum(eta + (Y - 1.) % Exp - Y % log(1. + phi) - lfac - frac);
  uword gh = w.size(), m = eta.size();
  vec phi = W * sigma, Exp(m), denom = 1. + phi;
  for(uword l = 0; l < gh; l++){
    vec this_eta = eta + tau * v[l];
    Exp += w[l] * log(trunc_exp(this_eta) + phi % Y);
  }
  return sum(
    (Y - 1.) % Exp- Y % log(1. + phi) - (trunc_exp(eta) + phi % Y)/denom
  );
}

// Update Gaussian by finding contribution from each individual i.d.
// Gaussian, variance = sigma (i.e. sigma \equiv simga^2_\epsilon)
// [[Rcpp::export]]
double sigma2_Gaussian_update(const arma::vec& eta, const arma::vec& Y, const arma::vec& tau, 
                              const arma::vec& w, const arma::vec& v){
  uword gh = w.size();
  double out = 0.0;
  for(uword l = 0; l < gh; l++){
    vec rhs = Y - eta - tau * v[l];
    out += w[l] * as_scalar(rhs.t() * rhs);
  }
  return out;
}

// The baseline hazard ----------------------------------------------------
// Update to \lambda (i.e. before usage in gamma update)
//' @keywords internal, this assumes mu_surv, tau_surv not calculated prior.
// [[Rcpp::export]]
arma::vec lambda_hat(Rcpp::List b, Rcpp::List Fu, Rcpp::List SS, Rcpp::List Sigma,
                     arma::vec& gamma_rep, arma::vec& zeta, 
                     arma::vec& nev, arma::vec& w, arma::vec& v){
  uword n = b.size(), gh = w.size(), unique_times = nev.size();
  vec out = zeros<arma::vec>(unique_times);
  mat g = diagmat(gamma_rep);
  
  for(uword i = 0; i < n; i++){
   
   vec b_i = b[i];
   mat SS_i = SS[i], Fu_i = Fu[i], S_i = Sigma[i];
   mat Q = Fu_i * g;
   mat A = Q * S_i * Q.t();
   
   // mu, tau.
   vec mu = exp(SS_i * zeta + Q * b_i), tau = sqrt(diagvec(A));
   
   int ui = Fu_i.n_rows;
   for(uword l = 0; l < gh; l++){
     out.subvec(0, ui - 1) += w[l] * mu % exp(tau * v[l]);
   }
   
  }
  return nev/out;
}

// Update to \lambda (i.e. after gamma updated)
//' @keywords internal
// [[Rcpp::export]]
arma::vec lambda_update(Rcpp::List b, Rcpp::List Fu, Rcpp::List SS, Rcpp::List Sigma,
                        Rcpp::List survtimes,
                        arma::vec& gamma_rep, arma::vec& zeta, 
                        arma::vec& nev, arma::vec& w, arma::vec& v){
  uword n = b.size(), gh = w.size(), unique_times = nev.size();
  vec out = zeros<arma::vec>(unique_times);
  mat g = diagmat(gamma_rep);
 
 for(uword i = 0; i < n; i++){
   
   vec b_i = b[i], st_i = survtimes[i];
   mat SS_i = SS[i], Fu_i = Fu[i], S_i = Sigma[i];
   mat Q = Fu_i * g;
   mat A = Q * S_i * Q.t();
   
   int num_survived = st_i.size();
   if(num_survived == 0) continue;
   
   // mu, tau.
   vec mu = exp(SS_i * zeta + Q * b_i), tau = sqrt(diagvec(A));
   
   int ui = Fu_i.n_rows;
   for(uword l = 0; l < gh; l++){
     out.subvec(0, ui - 1) += w[l] * mu % exp(tau * v[l]);
   }
   
 }
 return nev/out;
}

// (gamma, zeta) ----------------------------------------------------------
// [[Rcpp::export]]
double Egammazeta(const arma::vec& gammazeta, const arma::vec& b, const arma::mat& Sigma,
                  const arma::rowvec& S, const arma::mat& SS, const arma::mat& Fu, 
                  const arma::rowvec& Fi, const arma::vec& haz, const int Delta, 
                  const arma::vec& w, const arma::vec& v, const List& inds, const arma::uword K){
  uword gh = w.size(), q = Fu.n_cols;
  vec g = gammazeta.head(K), z = gammazeta.subvec(K, gammazeta.size() - 1), gammas(q);
  for(uword k = 0; k < K; k++){
    double gk = g[k];
    uvec inds_k = inds[k];
    gammas.elem(inds_k) += gk;
  }
  mat gmat = diagmat(gammas);
  mat Q = Fu * gmat;
  mat A = Q * Sigma * Q.t();
  vec mu = SS * z + Q * b, tau = sqrt(diagvec(A));
  double store = 0.;
  for(uword l = 0; l < gh; l++){
    store += w[l] * as_scalar(haz.t() * exp(mu + tau * v[l]));
  }
  return as_scalar(Delta * (S * z + Fi * (b % gammas))) - store;
}

//' @keywords internal
 // [[Rcpp::export]]
 arma::vec Sgammazeta(arma::vec& gammazeta, arma::vec& b, arma::mat& Sigma,
                      arma::rowvec& S, arma::mat& SS, arma::mat& Fu, arma::rowvec& Fi,
                      arma::vec& haz, int Delta, arma::vec& w, arma::vec& v, Rcpp::List b_inds, int K, double eps){
  int ps = gammazeta.size();
  vec out = vec(ps), epsvec = vec(ps);
  for(int i = 0; i < ps; i++){
   epsvec[i] = eps;
   double xi = 1.;//std::max(1.0, std::abs(gammazeta[i]));
   vec ge1 = gammazeta + xi * epsvec, ge2 = gammazeta - xi * epsvec;
   double f1 = Egammazeta(ge1, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K),
          f2 = Egammazeta(ge2, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K);
   out[i] = (f1 - f2) / (2 * eps);
   epsvec[i] = 0.;
  }
  return out;
}

//' @keywords internal
// [[Rcpp::export]]
arma::mat Hgammazeta(arma::vec& gammazeta, arma::vec& b, arma::mat& Sigma,
                     arma::rowvec& S, arma::mat& SS, arma::mat& Fu, arma::rowvec& Fi,
                     arma::vec& haz, int Delta, arma::vec& w, arma::vec& v, Rcpp::List b_inds, int K, double eps){
  int ps = gammazeta.size();
  mat out = zeros<mat>(ps, ps);
  vec epsvec = vec(ps, fill::value(eps));
  mat epsmat = diagmat(epsvec);
  double eps2 = pow(eps, 2.), eps42 = 4. * pow(eps, 2.);
  double f0 = Egammazeta(gammazeta, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K);
  for(int i = 0; i < (ps - 1); i++){
    vec eps_i = epsmat.col(i);
    // Diagonal terms
    vec ge_i1 = gammazeta + eps_i, ge_i2 = gammazeta - eps_i;
    double f1 = Egammazeta(ge_i1, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K),
           f2 = Egammazeta(ge_i2, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K);
    out(i,i) = (f1 - 2. * f0 + f2)/eps2;
    // Off-diagonal
    for(int j = (i + 1); j < ps; j++){
      vec eps_j = epsmat.col(j);
      vec ge_ij1 = gammazeta + eps_i + eps_j,
          ge_ij2 = gammazeta + eps_i - eps_j,
          ge_ij3 = gammazeta - eps_i + eps_j,
          ge_ij4 = gammazeta - eps_i - eps_j;
      double f1 = Egammazeta(ge_ij1, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K),
             f2 = Egammazeta(ge_ij2, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K),
             f3 = Egammazeta(ge_ij3, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K),
             f4 = Egammazeta(ge_ij4, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K);
      out(i,j) = (f1 - f2 - f3 + f4)/eps42;
      out(j,i) = out(i,j);
    }
  }
  // Calculate ps, psth item
  int last = ps-1;
  vec eps_i = epsmat.col(ps-1);
  vec ge_i1 = gammazeta + eps_i, ge_i2 = gammazeta - eps_i;
  double f1 = Egammazeta(ge_i1, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K),
         f2 = Egammazeta(ge_i2, b, Sigma, S, SS, Fu, Fi, haz, Delta, w, v, b_inds, K);
  out(last, last) = (f1 - 2. * f0 + f2)/eps2;
  return out;
}

// Metropolis algorithm for sampling from f(b_i|...; Omega) ---------------
//' @keywords internal
// [[Rcpp::export]]
List metropolis(const arma::vec& b, const List& Omega, 
                const List& Y, const List& X, const List& Z, const List& W,
                const List& family, const int Delta, const arma::rowvec& S, const arma::rowvec& Fi, 
                const double l0i, const arma::mat& SS, const arma::mat& Fu, const arma::rowvec& haz, 
                const arma::vec& gamma_rep, const List& beta_inds, const List& b_inds, const arma::uword K, 
                const arma::uword q, const int burnin, const int N, const arma::mat& Sigma, const double tune){
  // Unpack Omega
  mat D = Omega["D"];
  List sigma = Omega["sigma"];
  vec beta = Omega["beta"];
  vec zeta = Omega["zeta"];
  // MC stuff
  int iters = burnin + N;
  mat out = mat(q, iters);
  // Start
  int j = 1, num_accepts = 0;
  while(j < iters){
    double U = randu();
    vec b_current = out.col(j - 1);
    vec b_proposal = mvnrnd(b_current, tune * Sigma);
    double logf_current = -1. * joint_density(b_current, Y, X, Z, W, beta, D, sigma, family, Delta,
                                              S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double logf_proposal = -1. * joint_density(b_proposal, Y, X, Z, W, beta, D, sigma, family, Delta,
                                               S, Fi, l0i, SS, Fu, haz, gamma_rep, zeta, beta_inds, b_inds, K);
    double P = std::min(exp(logf_proposal - logf_current), 1.);
    if(U < P){
      b_current = b_proposal;
      if(j > burnin) num_accepts ++;
    }
    out.col(j) = b_current;
    j++ ;
  }
  
  // Remove burnin -----------
  out.shed_cols(0, burnin - 1);
  return List::create(_["walks"] = out,
                      _["burnin"] = burnin,
                      _["N"] = N,
                      _["AcceptanceRate"] = (double)num_accepts/(double)N);
}

// Fast dmvnorm and dmvt for use with dynpred-draws.R ---------------------
// (basically wrappers)
//' @keywords internal
// [[Rcpp::export]]
double dmvn_fast(const arma::vec& x, const arma::vec& mn,
                  const arma::mat& Sigma, const bool log__ = true){
   return as_scalar(
     dmvnrm_arma_fast(x.t(), mn.t(), Sigma, log__)
   );
 }

//' @keywords internal
// [[Rcpp::export]]
double dmvt_fast(const arma::vec& x, const arma::vec& mn,
                  const arma::mat& Sigma, const double df, const bool log__ = true){
   return as_scalar(
     dmvt_arma_fast(x.t(), mn.t(), Sigma, df, log__)
   );
 }
