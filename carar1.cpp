#include <TMB.hpp>

template <class Type>
Type f_sc(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_SCALAR(sd_phi);
  
  DATA_MATRIX(W);
  DATA_MATRIX(X);
  DATA_MATRIX(Y);    
  PARAMETER_VECTOR(beta);
  PARAMETER(sd);

  PARAMETER(inv_phi);
  //PARAMETER(sigma);
  
  
  Type phi_mean = 0.0;
  Type phi = f_sc(inv_phi);
  

  Type f=0;
  //f -= dnorm(phi, phi_mean, sd_pos, true);
  f -= dnorm(phi, phi_mean, sd_phi, true);

  

  Type variance = sd * sd;
  
  int n_t = Y.rows();
  int n_s = Y.cols();
  //std::cout << n_s <<"\n";
  
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
  
  
  vector<Type> Xb = X*beta;
  
  array<Type> residuals(n_t, n_s);
  
  for(int i=0; i < n_t; i++){
    for(int j=0; j<n_s; j++){
      residuals(i, j) = Y(i, j) - Xb(j);
    }
  }
  
  
  matrix<Type> I(n_s, n_s);
  matrix<Type> M(n_s, n_s);
  matrix<Type> C(n_s, n_s);
  C = W / n_s;
  I.setIdentity();
  M = variance * I;
  
  matrix<Type> Sigma(n_s, n_s);
  Sigma = (I - C).inverse() * M;
  


  f += SEPARABLE(MVNORM(Sigma), AR1(phi))(residuals);
  return(f);
  
}
