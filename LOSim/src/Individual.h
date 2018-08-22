/*
 * Individual.h
 *
 *  Created on: 31 Mar 2016
 *      Author: SH
 */

#ifndef INDIVIDUAL_H_INCLUDED__
#define INDIVIDUAL_H_INCLUDED__
#include <math.h>
//[[Rcpp::plugins(cpp11)]]

// [[RCPP::deppends(BH)]]
#include "Environment.h"
#include "RandomGenerator.h"
#include "Solar.h"
#include <boost/math/distributions.hpp>
#include <iostream>


// define pi which was dropped in C99
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

class Individual {
public:
	Individual(std::vector<double> currentParameters,
			std::vector<unsigned int> Observationtime,
			double startT,
			int startD,
			double solarLat, // latitude in Degree
			double solarLon, // longitude in Degree
 			double res); // resolution of grid cells
	virtual ~Individual();

	Solar* sol;
	double cellres;

	// hidden states
	std::vector<int> x;
	std::vector<int> y;
	std::vector<double> envVal;
	std::vector<double> timestamp;
	std::vector<double> stepDistance;
	std::vector<double> turningAngle;
	std::vector<double> energystamp;
	std::vector<int> activitystatus; // 0 = active dispersal, 1 = day-rest, 2 = 1 night + 2dayrests, ... n = n-1 nights + n dayrests
	std::vector<double> time;

	// observed states
	std::vector<double> xObs;
	std::vector<double> yObs;
	std::vector<double> envValObs;
	std::vector<double> rsc; // resource selection coefficient
	std::vector<double> stepDistanceObs;
	std::vector<double> turningAngleObs;
	std::vector<unsigned int> tObs; // in min // input variable

	bool observationModel;

	// parameters
	double habitatPreference;
	double stepShape;
	double stepScale;
	double directionalBias;
	double DispersalRestingTimeLogMean;
	double DispersalRestingTimeLogSd;
	double maximumEffort;
	double roostLambda;
	double observationError;
	int rscRange; // availability range in rsc_calc


	// intermediate states
	double energy;
	bool active;
	double daytime;
	double startDaytime;
	int doy;
	int startDoy;

	// validation checks
	unsigned int validityCode;
	unsigned int numberOfSteps;

	// movement model
	void move(unsigned int nSteps, int maximumDuration, Environment *env, RandomGenerator ranGen);
	// observation model
	void observe(RandomGenerator ranGen, Environment *env);
	// locate position where observed time matches simulated time
	unsigned int locate(unsigned int targetTime, unsigned int startPosition);
	// moore neighbourhood
	double mooreEnv(Environment *envir, int xcoord, int ycoord);
	// calculate resource selection coefficient
	double rsc_calc(Environment *envir, int xcoord, int ycoord);
	// update activity
	void updateActivity();
	// observe enery limit
	void checkEnergyLimit();
};


#endif /* INDIVIDUAL_H_ */
