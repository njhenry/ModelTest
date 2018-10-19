#include <TMB.hpp>

template <class Type>
Type f_sc(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(x);
  DATA_SCALAR(sd);

  PARAMETER(inv_phi);
  PARAMETER(sigma);
  
  
  Type phi_mean = 0.0;
  Type phi = f_sc(inv_phi);
  Type variance = sigma*sigma;
  

  Type f=0;
  //f -= dnorm(phi, phi_mean, sd_pos, true);
  f -= dnorm(phi, phi_mean, sd, true);
  f += SCALE(AR1(phi), variance)(x);
  //f += AR1(phi)(x);
  
  
  return(f);
}
