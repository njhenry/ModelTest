#include <RcppEigen.h>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>    

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

//function to calculate estimate of HSIC
// [[Rcpp::export]]
double fsic(NumericVector x, NumericVector y,
            NumericVector v, NumericVector w){
  //n is number of observations
  int n = x.size();
  int J = v.size();
  
  //allocate space for Gram matrices K, L and
  //centering matrix H
  Eigen::MatrixXd K(J, n);
  Eigen::MatrixXd L(J, n);
  
  //create Gram matrices
  //double scalex = mx * mx;
  for(int i = 0; i < J; i++){
    for(int j = 0; j < n; j++){
      //double val = exp(- pow (v[i] - x[j], 2) / scalex);
      double val = exp(- pow (v[i] - x[j], 2));
      K(i, j) = val;
    }
  }
  
  //double scaley = my * my;
  for(int i = 0; i < J; i++){
    for(int j = 0; j < n; j++){
      //double val = exp(- pow (w[i] - y[j], 2) / scaley);
      double val = exp(- pow (w[i] - y[j], 2));
      L(i, j) = val;
    }
  }
  
  Eigen::VectorXd u;
  
  u = (1.0/(n-1))*(K.cwiseProduct(L)).rowwise().sum() - (1.0/(n*(n-1)))*K.rowwise().sum().cwiseProduct(L.rowwise().sum());
  
  Eigen::VectorXd ub;
  Eigen::MatrixXd G(J, n);
  
  ub = (1.0/n)*K.cwiseProduct(L).rowwise().sum() - (1.0/pow(n, 2))*K.rowwise().sum().cwiseProduct(L.rowwise().sum());
  G = (K - (1.0/n)*K.rowwise().sum().replicate(1, n)).cwiseProduct(L - (1.0/n)*L.rowwise().sum().replicate(1, n)) -
    ub.replicate(1, n);
  //  L - (1.0/n)*(L.rowwise().sum()).replicate(1, n));
  double gn = 0.00001;
  Eigen::MatrixXd Sdash  = (1.0/n)*G * G.transpose() + gn * MatrixXd::Identity(J, J);
  
  Eigen::VectorXd temp = Sdash.householderQr().solve(u);
  
  double ln = n * u.dot(temp);
  return ln;
}


//function to do bootstrap and get p-value for FSIC
// [[Rcpp::export]]
double fstest(NumericVector x, NumericVector y,
              NumericVector v, NumericVector w, int N_samples){
  srand(time(0));
  double test_stat = fsic(x, y, v, w);
  
  NumericVector x2 = clone(x);
  std::vector<double> dist(N_samples);
  
  
  for(int i = 0; i < N_samples; i++){
    std::random_shuffle(x2.begin(), x2.end());
    dist[i] = fsic(x2, y, v, w);
    std::cout << dist[i] <<"\n";
  }
  std::sort(dist.begin(), dist.end());
  
  if(test_stat <= dist[0]){
    return(1);
  }
  
  int i = 1;
  while((test_stat > (dist[i] + 0.001)) && (i < N_samples)){
    std::cout << "\n" << i << "\n";
    std::cout << dist[i] << "  " << test_stat << "\n";
    std::cout << N_samples << "\n";
    i += 1;
    std::cout << (double) (test_stat >= dist[i]) << "\n";
  }
  std::cout << "\n" << dist[i] << "\n";
  std::cout <<"\n"<<test_stat<<"\n";
  std::cout << i << "\n";
  return(1 - (double) (i - 1) / (double) N_samples);
}

// [[Rcpp::export]]
double fstest_debug(NumericVector x, NumericVector y,
              NumericVector v, NumericVector w, int N_samples){

  double test_stat = fsic(x, y, v, w);
  
  NumericVector x2 = clone(x);
  std::vector<double> dist(N_samples);
  
  
  for(int i = 0; i < N_samples; i++){
    srand(time(0));
    std::random_shuffle(x2.begin(), x2.end());
    dist[i] = fsic(x2, y, v, w);
    std::cout << dist[i] <<"\n";
  }
  std::sort(dist.begin(), dist.end());
  
  if(test_stat <= dist[0]){
    return(1);
  }
  
  int i = 1;
  while((test_stat > dist[i]) && (i < N_samples)){
    std::cout << "\n" << i << "\n";
    std::cout << dist[i] << "  " << test_stat << "\n";
    std::cout << N_samples << "\n";
    i += 1;
    std::cout << (double) (test_stat >= dist[i]) << "\n";
  }
  std::cout << "\n" << dist[i] << "\n";
  std::cout <<"\n"<<test_stat<<"\n";
  std::cout << i << "\n";
  return(1 - (double) (i - 1) / (double) N_samples);
}
