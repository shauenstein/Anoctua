#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector rsc_calc_obs(NumericMatrix environment, 
                            IntegerMatrix cellcoords,
                            int rscRange) {
  
  IntegerVector rscLegs(2 * rscRange + 1); // array of opposite legs
  int rscIndex = 0; // loop subscript
  for(int aR = - rscRange; aR <= rscRange; aR++, rscIndex++){
    // get indices for area within perception
    rscLegs(rscIndex) = floor(sqrt(std::pow(rscRange, 2.0) - std::pow(aR, 2.0)));
  }
  
  NumericVector out(cellcoords.nrow());
  for(int i = 0; i < cellcoords.nrow(); i++){
    if(i == 0){
      out(i) = NumericVector::get_na();
    }else{
      double minEnv = 1;
      double maxEnv = 0;
      int sEnv = 1; // number of env. values smaller than at current location
      int lEnv = 1; // number of env. values larger than at current location
      double curEnv = environment(cellcoords(i,1), cellcoords(i,0)); // current env. value
  
      // reset subscript
      rscIndex = 0;
      // loop through cells within rscRange
      for(int h = - rscRange; h <= rscRange; h++, rscIndex++){
        // vertical loop, y-direction
        for(int v = rscLegs[rscIndex]; v >= -rscLegs[rscIndex]; v--){
  
          // skip point of origin, give a weight of zero
          if(h != 0 && v != 0){
            // check if env is larger or smaller than curEnv
            double testEnv = environment(cellcoords(i,1)+v, cellcoords(i,0)+h);
            if(testEnv > curEnv){
              lEnv += 1;
            }else{
              if(testEnv < curEnv){
                sEnv += 1;
              }
            }
            // update min and max env. values
            if(testEnv < minEnv){
              minEnv = testEnv;
            }
            if(testEnv > maxEnv){
              maxEnv = testEnv;
            }
          }
        }
      }
  
      // ratio divided by range
      if((maxEnv - minEnv) == 0){
        out(i) = NumericVector::get_na();
      }else{
        // out(i) = (lEnv/(double) sEnv) / (maxEnv - minEnv);
        out(i) = lEnv/(double) sEnv;
      }
    }
    
  }
  
  return out;
}


