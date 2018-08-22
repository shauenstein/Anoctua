#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector dObsSim(NumericMatrix simulated, 
                            NumericVector observed,
                            NumericVector reference) {
  int n = simulated.nrow();
  NumericVector distances(n);
  
  for(int i = 0;i < n; i++){
    distances(i) = sqrt(sum(pow((simulated(i,_) - observed) / reference,2)));
  }
  return distances;
}


