#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(response)
  DATA_STRUCT(spde,spde_t);
  DATA_STRUCT(spde_2,spde_t);
  DATA_VECTOR(cov);
  DATA_VECTOR(cov_2);
  DATA_SPARSE_MATRIX(A);
  
  PARAMETER_VECTOR(gp);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(gp_2);
  PARAMETER(log_kappa_2);
  PARAMETER(log_tau_2);
  PARAMETER(sigma);
  
  
  
  Type kappa = exp(log_kappa);
  Type kappa_2 = exp(log_kappa_2);
  Type sigma_mean = 1.0;
  Type sigma_sd = 1.0;
  Type log_kappa_mean = 1.0;
  Type log_kappa_sd = 1.0;
  Type log_tau_mean = 0.0;
  Type log_tau_sd = 1.0;
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  SparseMatrix<Type> Q_2 = Q_spde(spde, kappa_2);
  int n = response.size();
  
  Type f;
  f = 0;
  f += SCALE(GMRF(Q), 1/exp(log_tau))(gp);
  f += SCALE(GMRF(Q_2), 1/exp(log_tau_2))(gp_2);
  f -= dnorm(log_kappa, log_kappa_mean, log_kappa_sd, true);
  f -= dnorm(log_tau, log_tau_mean, log_tau_sd, true);
  f -= dnorm(log_kappa_2, log_kappa_mean, log_kappa_sd, true);
  f -= dnorm(log_tau_2, log_tau_mean, log_tau_sd, true);
  f -= dnorm(sigma, sigma_mean, sigma_sd, true);
  
  vector<Type> field = A * gp;
  vector<Type> field_2 = A * gp_2;
  vector<Type> responsepred(n);
  
  for(int i=0; i<n; i++){
    responsepred[i] = field[i] * cov[i] + field_2[i] * cov_2[i];
  }
  
  for(int i=0; i<n; i++){
    f -= dnorm(responsepred[i], response[i], sigma, true);
  }
  
  return(f);
}
