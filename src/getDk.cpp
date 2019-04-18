#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <cmath>
#include "include/getDk.hpp"
using namespace Rcpp;
using namespace arma;

//' Discrete derivative matrix
//'
//' \code{get_D1} computes the sparse discrete derivative matrix.
//'
//' @param n length of input
//' @export
// [[Rcpp::export]]
arma::sp_mat get_D1(int n){
  int numberNonZero = 2*(n-1);
  arma::vec values = ones<vec>(numberNonZero);
  values.subvec(n-1, 2*(n-1)-1) = -1*values.subvec(n-1, 2*(n-1)-1);
  arma::umat locs = repmat(linspace<urowvec>(0,n-2,n-1),2,2);
  locs.submat(1, n-1, 1, numberNonZero-1) = locs.submat(1, n-1, 1,
              numberNonZero-1) + 1;
  arma::sp_mat D1 = arma::sp_mat(locs, values);
  return D1;
}

//' kth order sparse difference matrix
//'
//' \code{get_Dkn} computes the sparse discrete kth derivative matrix
//'
//' @param n length of input
//' @param k order of the derivative
//' @export
// [[Rcpp::export]]
arma::sp_mat get_Dk(int n,
                    int k){
  arma::sp_mat D = get_D1(n);
  for (int i=2; i < k+1; i++){
    D = get_D1(n-i+1)*D;
  }
  // Rcout << "D_k" << std::endl << D << std::endl;
  return D;
}
