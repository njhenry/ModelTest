#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(Y);
  DATA_MATRIX(X);
  DATA_MATRIX(W);
  DATA_SCALAR(phi_min);
  DATA_SCALAR(phi_max);
  //DATA_MATRIX(In);
  PARAMETER(sigma);
  PARAMETER(phi);
  PARAMETER_VECTOR(beta);
  
  int n = Y.size();
  matrix<Type> C = phi * W;
  Type f = 0;
  
  /*
  //priors
  //gamma prior for precision
  Type shape = 3.0;
  Type scale = 1.0;
  f -= dgamma(1.0/sigma, shape, scale, true);
  //rescale phi to be between zero and one
  Type phi_scaled = (phi - phi_min) / (phi_max - phi_min);
  Type phi_alpha = 2.0;
  Type phi_beta = 2.0;
  f -= dbeta(phi_scaled, phi_alpha, phi_beta, true);
  */
  
  matrix<Type> I(n,n);
  I.setIdentity();
  
  matrix<Type> S;
  S = (I - C).inverse() * (sigma * I);
  f += MVNORM(S)(Y - X*beta);
  
  return(f);
}
