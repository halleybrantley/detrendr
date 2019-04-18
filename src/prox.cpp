#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include "include/getDk.hpp"
#include "include/prox.hpp"
using namespace Rcpp;
using namespace arma;

//' @useDynLib detrendr
//' @importFrom Rcpp evalCpp

//' Banded Cholesky Solve
//'
//' \code{chol_solve} Solves a linear system cholM%*%x=b
//' when cholM is a sparse banded cholesky.
//'
//' @param cholM sparse banded cholseky decomposition of discrete derivative
//' matrix of order k
//' @param b dense solution vector
//' @param k order of discrete derivative matrix
//' @param upper boolean indicator of whether cholM is upper or lower triangular
//' @export
// [[Rcpp::export]]
arma::vec chol_solve(arma::sp_mat cholM,
                     arma::vec b,
                     int k,
                     bool upper = true){
  int n = cholM.n_cols;
  arma::uvec ind;
  arma::vec x=b;

  if(upper){
    ind = regspace<uvec>(n-1, 0);
  } else {
    ind = regspace<uvec>(0, n-1);
  }

  for (int i=0; i < k+1; i++){
    for (int j=0; j < i; j++){
      x(ind(i)) -= cholM(ind(i), ind(j))*x(ind(j));
    }
    x(ind(i)) /= cholM(ind(i), ind(i));
  }

  for (int i=k+1; i < n; i++){
    for (int j = i-k; j < i; j++){
      x(ind(i)) -= cholM(ind(i), ind(j))*x(ind(j));
    }
    x(ind(i)) /= cholM(ind(i), ind(i));
  }
  return x;
}


//' Proximal Mapping
//'
//' \code{prox_quantile} computes the proximal mapping of the check function.
//'
//' @param w input
//' @param tau quantile parameter
//' @param alpha scale parameter
//' @export
// [[Rcpp::export]]
arma::vec prox_quantile(arma::vec w,
                        double tau,
                        double alpha){
  int n = w.n_elem;
  arma::vec prox_out = zeros<vec>(n);
  double threshold1 = tau*alpha;
  double threshold2 = -(1 - tau)*alpha;

  for (int i=0; i<n; i++){
    if (w(i) > threshold1) {
      prox_out(i) = w(i) - threshold1;
    } else if (w(i) < threshold2){
      prox_out(i) = w(i) - threshold2;
    }
  }
  return prox_out;
}



//' Proximal mapping of f_1
//'
//' \code{prox_f1} computes the proximal mapping of the average quantile loss
//'
//' @param theta input
//' @param y response
//' @param tau quantile parameter
//' @param step step-size
//' @export
// [[Rcpp::export]]
arma::vec prox_f1(arma::vec theta,
                  arma::vec y,
                  double tau = 0.05,
                  double step = 1.0){
  //int n = theta.n_elem;
  arma::vec w = y - theta;
  return y - prox_quantile(w, tau, step);
}


//' Proximal mapping of f_2
//'
//' \code{prox_f2} computes the proximal mapping of the L1 penalty
//'
//' @param eta input
//' @param lambda regularization parameter
//' @param step step-size
//' @examples
//' set.seed(12345)
//' n <- 1e3
//' eta <- seq(-3, 3, length.out=n)
//' lambda <- 1
//' prox_out <- prox_f2(eta, lambda)
//' plot(eta, prox_out, type = 'l')
//' abline(0,1)
//' @export
// [[Rcpp::export]]
arma::vec prox_f2(arma::vec eta,
                  double lambda,
                  double step = 1){
  return prox_quantile(eta, 0.5, 2*step*lambda);
}


//' Proximal mapping
//'
//' \code{prox} computes the block separable proximal mapping, changes theta
//' and eta in place
//'
//' @param theta input
//' @param eta input
//' @param y response
//' @param lambda regularization parameter
//' @param tau quantile parameter
//' @param step step-size
//' @export
// [[Rcpp::export]]
void prox(arma::vec& theta,
                arma::vec& eta,
                arma::vec y,
                double lambda,
                double tau = 0.05,
                double step = 1.0){
  theta = prox_f1(theta, y, tau, step);
  eta = prox_f2(eta, lambda, step);
}



//' Project onto subspace
//'
//' \code{project_V} projects (theta, eta) onto the subspace eta = D*theta.
//' Updates values of theta and eta in place.
//'
//' @param theta first input
//' @param eta second input
//' @param D differencing matrix
//' @param cholM upper triangular cholesky decomposition of  I + DtD
//' @param k order of differencing matrix
//' @export
// [[Rcpp::export]]
void project_V(arma::vec& theta,
                     arma::vec& eta,
                     arma::sp_mat D,
                     arma::sp_mat cholM,
                     int k){
  arma::vec DtEta = vectorise(D.t()*eta);
  theta = chol_solve(cholM, chol_solve(cholM.t(), theta + DtEta, k, false),
                     k, true);
  eta = D*theta;
}


