#ifndef __prox__
#define __prox__


double check_loss(arma::vec r,
                  double tau);

arma::vec chol_solve(arma::sp_mat cholM,
                     arma::vec b,
                     int k,
                     bool upper);

arma::vec prox_quantile(arma::vec w,
                        double tau,
                        double alpha);

arma::vec prox_f1(arma::vec theta,
                  arma::vec y,
                  double tau,
                  double step);

arma::vec prox_f2(arma::vec eta,
                  double lambda,
                  double step);

void prox(arma::vec& theta,
          arma::vec& eta,
          arma::vec y,
          double lambda,
          double tau,
          double step);
  
void project_V(arma::vec& theta,
               arma::vec& eta,
               arma::sp_mat D,
               arma::sp_mat cholM,
               int k);  

#endif // __prox__