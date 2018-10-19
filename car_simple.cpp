#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_MATRIX(W);
  DATA_MATRIX(X);
  DATA_VECTOR(Y);
  PARAMETER_VECTOR(beta);
  PARAMETER(sd);
  
  int n = Y.size();
  matrix<Type> C = W / n;
  vector<Type> Xb = X*beta;
  vector<Type> residuals = Y - Xb;
  vector<Type> res_sum = C * residuals;
  Type variance = sd * sd;
  Type sd_pos = sqrt(variance);
  
  Type f=0;
  
  //priors
  Type beta_mean = 0.0;
  Type beta_sd = 1.0;
  Type precision = 1.0 / variance;
  Type pr_sh = 3.0;
  Type pr_sc = 3.0;
  
  
  
  
  int n_beta = beta.size();
  for(int i=0; i<n_beta; i++){
    f -= dnorm(beta(i), beta_mean, beta_sd, true);
  }
  f -= dgamma(precision, pr_sh, pr_sc, true);
   
  
  
  for(int i=0; i<n; i++){
    //f -= dnorm(Y(i), Xb(i) + res_sum(i), variance, true);
    f -= dnorm(Y(i), Xb(i), sd_pos, true);
    //f -= dnorm(Y(i), Xb(i), sd);
  }

  
  return(f);
}
