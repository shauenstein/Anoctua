//============================================================================
// Name        : ST2.cpp
// Author      : SH
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>

#include "Model.h"
#include "Environment.h"
#include "RandomGenerator.h"
#include "Solar.h"

// define pi which was dropped in C99
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif


// function to print matrix as is
void displayMatrix(std::vector<std::vector<unsigned int> > dArray, unsigned int xSize,
		unsigned int ySize){

	std::cout << "The Matrix looks like this:\n" << std::endl;

//	std::vector<double>* pArray = dArray;
	unsigned int subscript = 1; // running subscript
	for(unsigned int y = 0; y < ySize; y++){
		for(unsigned int x = 0; x < xSize; x++, subscript++){
			std::cout << dArray[x][y] << " ";
			if(subscript % xSize == 0){
				std::cout << std::endl;
			}
		}
	}
	std::cout << std::endl;
}

// equivalent to dnorm() in R
double gaussianPdf(double x, double mean = 0.0, double sd = 1.0, bool Log = true){
	double pOut;
	pOut = 1 / (sd * sqrt(2 * M_PI)) * exp(-pow(x - mean, 2) / (2 * pow(sd, 2)));
	if(Log){
		pOut = log(pOut);
	}
	return pOut;
}


// main // to be replaced with Rcpp chunk
int main() {
//	Solar sol;
//	sol.solarcalc(31, 47.999, 7.8421);
//	std::cout <<  "sunrise = " << sol.sunrise <<
//			" sunset = " << sol.sunset << std::endl;

	// read txt matrix
//    std::ifstream in("D:/InteamDocs/Uni/littleowl/test.txt", std::ios_base::in | std::ios_base::binary);
//    std::vector<std::vector<double> >*v = new std::vector<std::vector<double> >
//    	 (0, std::vector<double> (0));
//    if (in) {
//        std::string line;
//
//        while (std::getline(in, line)) {
//          std::istringstream buffer(line);
//          std::vector<double> line((std::istream_iterator<double>(buffer)),
//                                   std::istream_iterator<double>());
//
//          v->push_back(line);
//        }
////        while (std::getline(in, line)) {
////            v->push_back(std::vector<double>());
////
////            // Break down the row into column values
////            std::stringstream split(line);
////            double value;
////
////            while (split >> value)
////                v->back().push_back(value);
////        }
//    }
//
//    for (int i = 0; i < v->size(); i++) {
//        for (int j = 0; j < v->at(i).size(); j++)
//            std::cout << v->at(i).at(j) << ' ';
//
//        std::cout << '\n';
//    }



	//
	RandomGenerator ranGen;

	// use std::vector
	unsigned int lat = 10000; // latitude, height
	unsigned int lon = 10000; // longitude, width
	unsigned int randomSeed = 15;

	std::vector<double> param;
	unsigned int startX = 5000;
	unsigned int startY = 5000;
	param.push_back(startX);
	param.push_back(startY);
	param.push_back(1); // habitat preference
	param.push_back(2.8); // stepShape
	param.push_back(2); // stepScale
	param.push_back(1.5); // directional bias
	param.push_back(1.5); // dispersal resting time logMean
	param.push_back(0.7); // dispersal resting time logSd
	param.push_back(10); // maximum effort
	param.push_back(2); // roost lambda
//	param.push_back(4); // recoveryRate
//	param.push_back(-5); // beta1
//	param.push_back(5); // beta2
	param.push_back(1); // observationError
	param.push_back(5); // rsc range


	double startT = 14.32;
	int startDoy = 150;
	double solLat = 47.999;
	double solLon = 7.8421;
	double resolution = 1;

	std::vector<unsigned int> times;
	times.push_back(0);
	times.push_back(2373);
	times.push_back(6863);
	times.push_back(6868);
	times.push_back(13909);
	times.push_back(13914);
	times.push_back(21381);
	times.push_back(21386);
	times.push_back(22030);
	times.push_back(27039);
	times.push_back(27044);
	times.push_back(34218);
	// = {0, 2373, 2378, 6863,
	// 		6868, 13909, 13914, 21381, 21386, 22030, 27039,
	// 		27044, 34218};


	std::vector<std::vector<double> >*landscape = new std::vector<std::vector<double> >
	 (lat, std::vector<double> (lon));

	// fill array, row major
//	for(unsigned int y = 0; y < lon; y++){
//		for(unsigned int x = 0; x < lat; x++){
//			landscape->at(x).at(y) = (gaussianPdf(y, lon/2, 2, true) +
//										 gaussianPdf(x, lat/2, 2, true)) * -1;
//		}
//	}

	for(unsigned int y = 0; y < lon; y++){
		for(unsigned int x = 0; x < lat; x++){
			landscape->at(x).at(y) = 1;
		}
	}

	// set random seed
	ranGen.seedrand(randomSeed);

//	bool checkBinom = ranGen.rBinomial(0.9);
//	std::cout << checkBinom << std::endl;

	Model mod(landscape);
	//displayMatrix(*mod.env->values, lat,lon);
	std::vector<unsigned int> tt;
	int maxDuration = -999;
	Individual* indOut = mod.runSimulation(param, 2000, maxDuration, ranGen, times,
			startT, startDoy, solLat, solLon, resolution);
//	std::cout << indOut->observationIndex.at(0)
//			<< "\n" << indOut->observationIndex.at(1)
//			<< "\n" << indOut->observationIndex.at(2)
//			<< "\n" << indOut->observationIndex.at(3)
//			<< " and nSim: " << indOut->x.size()
//			<< std::endl;
//	std::cout << indOut->x.at(1) << " and rsc2: " << indOut->rsc.at(2) << std::endl;
//	std::cout << indOut->validityCode << " " << indOut->numberOfSteps << std::endl;
	//displayMatrix(out, 2,6);
	delete landscape;

	// compute 95% quantile for gamma distribution
//	double shape = 2;
//	double scale = 1;
//	boost::math::gamma_distribution<double> d(shape, scale);
//	std::cout << quantile(d, 0.95) << std::endl;


	return 0;
}

