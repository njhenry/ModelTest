#include <TMB.hpp>



template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_MATRIX(W);
  DATA_MATRIX(X);
  DATA_VECTOR(Y);    
  PARAMETER_VECTOR(beta);
  PARAMETER(phi);
  PARAMETER(sd);
  Type variance = sd * sd;
  
  int n = Y.size();
  Type f = 0;
  
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
  
  matrix<Type> I(n, n);
  matrix<Type> M(n, n);
  matrix<Type> C(n, n);
  C = phi * W;
  I.setIdentity();
  M = variance * I;
  
  matrix<Type> Sigma;
  Sigma = (I - C).inverse() * M;

  
  f += MVNORM(Sigma)(Y - X*beta);
  return(f);
}
