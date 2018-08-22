
#include <Rcpp.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


#include "Environment.h"
#include "Individual.h"
#include "Model.h"
#include "RandomGenerator.h"
#include "Solar.h"

using namespace Rcpp;

// [[Rcpp::export]]
List callModel(NumericMatrix environment, unsigned int iterations, 
               int maxDuration, NumericMatrix parameters, 
               unsigned int randomSeed, int output, 
               IntegerVector TimeAtObservation, NumericVector upperleft, 
               double environmentResolution, Rcpp::CharacterVector filename,
               double startTime, int startDayOfYear, 
               double solarLatitude, double solarLongitude) {
  
  
  // initialise random generator and random seed
  RandomGenerator ranGen;
  ranGen.seedrand(randomSeed);

  // convert environment Rcpp::NumericMatrix to std::vector<std::vector>>
  // or read from file
  int nr = environment.nrow(), nc = environment.ncol();
  std::vector<std::vector<double> >* landscape = new std::vector<std::vector<double> >
    (nr, std::vector<double> (nc));
  
  if(filename[0] == "empty"){
    for(int c = 0; c < nc; c++){
      for(int r = 0; r < nr; r++){
        landscape->at(r).at(c) = environment(r, c);
      }
    }
  }else{
    std::string fname = as<std::string>(filename[0]);
    // read txt matrix
    std::ifstream inFile(fname);
    
    if (inFile.good()) {
      std::string temp;
      
      while (std::getline(inFile, temp)) {
        std::istringstream buffer(temp);
        std::vector<double> line((std::istream_iterator<double>(buffer)),
                                 std::istream_iterator<double>());
        
        landscape->push_back(line);
      }
    }
    
    inFile.close();
  }
  
  // create output lists
  Rcpp::List simulation = Rcpp::List::create();
  
  // validity vector
  Rcpp::IntegerVector validity(parameters.nrow());
  // end vector
  Rcpp::IntegerVector nSteps(parameters.nrow());
  
  // create model
  Model mod(landscape);
  
  
  // initiate current observation time series, if !ABCcalibration it stays empty
  std::vector<unsigned int> observationTimeVec;
  // Do we calibrate?
  if(output != 1){
    // transform current row in TimeAtObservation to std::vector
    for(int t = 0; t < TimeAtObservation.size(); t++){
      observationTimeVec.push_back((unsigned int) TimeAtObservation(t));
    }
  }
  
  for(int p = 0; p < parameters.nrow(); p++){

    // transform current row in parameters to std::vector
    std::vector<double> currentParameters;
    for(int t = 0; t < parameters.ncol(); t++){
      currentParameters.push_back(parameters(p,t));
    }
    // std::cout << "C1: x0 = " << currentParameters.at(0) << "; y0 = " << currentParameters.at(1) << std::endl;
    Individual* indOut = mod.runSimulation(currentParameters, iterations, maxDuration, ranGen, observationTimeVec,
                                           startTime, startDayOfYear, solarLatitude,solarLongitude, environmentResolution); 
    // create dataframe fitting the output requirements

    // return all simulated steps without observational imprecision
    if(output != 2){
      NumericMatrix outMat (indOut->x.size(), 3);
      for(int i = 0; i < outMat.nrow(); i++){
        outMat(i, 0) = indOut->timestamp.at(i);
        outMat(i, 1) = indOut->x.at(i) * environmentResolution + upperleft(0);
        outMat(i, 2) = -indOut->y.at(i) * environmentResolution + upperleft(1);
        // outMat(i, 3) = indOut->stepDistance.at(i) * environmentResolution;
        // if(i < 2){
        //   outMat(i, 4) = NA_REAL;
        // }else{
        //   outMat(i, 4) = indOut->turningAngle.at(i);
        // }
        // outMat(i, 5) = indOut->envVal.at(i);
        // outMat(i, 6) = indOut->energystamp.at(i);
        // outMat(i, 7) = indOut->activitystatus.at(i);
        // outMat(i, 8) = indOut->time.at(i);
        
      }
      
      colnames(outMat) = CharacterVector::create("timestamp", "x", "y");//, "stepDistance", 
               //"turningAngle", "habitatSuitability", "energylevel", "activitystatus", "daytime");
      
      simulation.push_back(outMat);
      
    }else{
      // return observed
      if(output != 1){
        NumericMatrix outMat (indOut->xObs.size(), 7);
        for(int i = 0; i < outMat.nrow(); i++){
          outMat(i, 0) = observationTimeVec.at(i);
          outMat(i, 1) = indOut->xObs.at(i) * environmentResolution + upperleft(0);
          outMat(i, 2) = -indOut->yObs.at(i) * environmentResolution + upperleft(1);
          outMat(i, 3) = indOut->stepDistanceObs.at(i) * environmentResolution;
          if(i < 2){
            outMat(i, 4) = NA_REAL;
          }else{
            outMat(i, 4) = indOut->turningAngleObs.at(i);
          }
          outMat(i, 5) = indOut->envValObs.at(i);
          outMat(i, 6) = indOut->rsc.at(i); 
        }
        
        colnames(outMat) = CharacterVector::create("timestamp", "xObserved", "yObserved", 
                 "stepDistanceObserved", "turningAngleObserved", //"habitatSuitability",
                 "habitatSuitabilityObserved", "rsc"); //,"energylevel", "activitystatus", "daytime");
        
        simulation.push_back(outMat);
      
      // raw output with boolean indication where observed would be
      }else{
        // NumericMatrix outMat (indOut->x.size(), 12);
        // for(int i = 0; i < outMat.nrow(); i++){
        //   outMat(i, 0) = indOut->timestamp.at(i);
        //   outMat(i, 1) = indOut->x.at(i) * environmentResolution + upperleft(0);
        //   outMat(i, 2) = -indOut->y.at(i) * environmentResolution + upperleft(1);
        //   if(std::isnan(indOut->xObs.at(i))){
        //     outMat(i, 3) = NA_REAL;
        //     outMat(i, 4) = NA_REAL;
        //     outMat(i, 8) = NA_REAL;
        //     outMat(i, 9) = 0;
        //   }else{
        //     outMat(i, 3) = indOut->xObs.at(i) * environmentResolution + upperleft(0);
        //     outMat(i, 4) = -indOut->yObs.at(i) * environmentResolution + upperleft(1);
        //     outMat(i, 8) = indOut->envValObs.at(i);
        //     outMat(i, 9) = 1;
        //   }
        //   if(i == 0){
        //     outMat(i, 5) = NA_REAL;
        //   }else{
        //     outMat(i, 5) = indOut->stepDistance.at(i) * environmentResolution;
        //   }
        //   if(i < 2){
        //     outMat(i, 6) = NA_REAL;
        //   }else{
        //     outMat(i, 6) = indOut->turningAngle.at(i);
        //   }
        //   outMat(i, 7) = indOut->envVal.at(i);
        //   outMat(i, 10) = indOut->energystamp.at(i);
        //   outMat(i, 11) = indOut->activitystatus.at(i);
        // }
        // 
        // 
        // colnames(outMat) = CharacterVector::create("timestamp", "x", "y", "xObserved", "yObserved",
        //          "stepDistance", "rel.turningAngle", "habitatSuitability", "habitatSuitabilityObserved",
        //          "observed", "energylevel", "activitystatus");
        // 
        // simulation.push_back(outMat);
        // 
      }
    }
    validity(p) = indOut->validityCode;
    nSteps(p) = indOut->numberOfSteps;
    delete indOut;

  // end parameter combination  
  }
  
  delete landscape;
  
  return Rcpp::List::create(_["simulation"] = simulation,
                            _["validity"] = validity,
                            _["nSteps"] = nSteps);
                              
}

