#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include "include/getDk.hpp"
#include "include/prox.hpp"
using namespace Rcpp;
using namespace arma;

//'
//' One step of Spingarn's algorithm
//'
//' \code{spingarn_one_step} updates theta and eta in place
//' @param theta input 1
//' @param eta input 2
//' @param y response
//' @param D differencing matrix
//' @param cholM upper cholesky of  (I + DtD)
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @param k order of differencing matrix
//' @export
// [[Rcpp::export]]
void spingarn_one_step(arma::vec& theta,
                       arma::vec& eta,
                       arma::vec y,
                       arma::sp_mat D,
                       arma::sp_mat cholM,
                       double lambda,
                       double tau = 0.05,
                       double step = 1,
                       int k = 3){
  arma::vec theta_old = theta;
  arma::vec eta_old = eta;
  prox(theta, eta, y, lambda, tau, step);
  arma::vec thetaMid = 2*theta-theta_old;
  arma::vec etaMid = 2*eta-eta_old;
  project_V(thetaMid, etaMid, D, cholM, k);
  theta = theta_old + 1*(thetaMid - theta);
  eta = eta_old + 1*(etaMid - eta);
}

//'
//' Multiple steps of Spingarn's algorithm
//'
//' \code{spingarn_multi_step}
//' @param theta input 1
//' @param eta input 2
//' @param y response
//' @param D discrete differencing matrix
//' @param cholM cholesky of I + DtD
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @param numberIter number of iterations
//' @param k order of differencing
//' @export
// [[Rcpp::export]]
Rcpp::List spingarn_multi_step(arma::vec theta,
                               arma::vec eta,
                               arma::vec y,
                               arma::sp_mat D,
                               arma::sp_mat cholM,
                               double lambda,
                               double tau = 0.05,
                               double step = 1,
                               double numberIter=1,
                               int k=3){
  arma::vec Vdiff = zeros<vec>(1);
  arma::vec theta_cp = theta;
  arma::vec eta_cp = eta;
  
  for (int i = 1; i < numberIter; i++){
    spingarn_one_step(theta_cp, eta_cp, y, D, cholM,
                      lambda, tau, step, k);
    if (i % 100 == 0){
      Rcpp::checkUserInterrupt();
    }
  }
  
  prox(theta_cp, eta_cp, y, lambda, tau, step);
  return Rcpp::List::create(Named("theta")=theta_cp,
                            _["eta"]=eta_cp);
}


