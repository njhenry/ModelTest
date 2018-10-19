#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_MATRIX(W);
  DATA_MATRIX(X);
  DATA_VECTOR(Y);
  PARAMETER_VECTOR(beta);
  PARAMETER(variance);
  
  int n = Y.size();
  matrix<Type> C = W / n;
  std::cout << C.sum()<<"\n";
  vector<Type> Xb = X*beta;
  vector<Type> residuals = Y - Xb;
  vector<Type> res_sum = C * residuals;
  
  Type f=0;
  
  //priors
  Type beta_mean = 0.0;
  Type beta_sd = 1.0;
  Type precision = 1.0 / variance;
  Type pr_sh = 1.0;
  Type pr_sc = 0.2;
  
  
  int n_beta = beta.size();
  for(int i=0; i<n_beta; i++){
    f -= dnorm(beta(i), beta_mean, beta_sd, true);
  }
  f -= dgamma(precision, pr_sh, pr_sc, true);
  
  
  for(int i=0; i<n; i++){
    f -= dnorm(Y(i), Xb(i) + res_sum(i), variance, true);
  }
  std::cout << Xb(0) << "\n";
  
  REPORT(f);
  return(f);
}
