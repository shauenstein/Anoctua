/*
 * Model.cpp
 *
 *  Created on: 14 Mar 2016
 *      Author: SH
 */

#include "Model.h"



Model::Model(
		// environments input
		std::vector<std::vector<double> >*habitatSuitability) {
	env = new Environment(habitatSuitability);

	// TODO Auto-generated constructor stub
}

Model::~Model() {
	// TODO Auto-generated destructor stub
}



// run simulation
Individual* Model::runSimulation(std::vector<double> currentParameters,
		unsigned int iterations, int maxDuration,
		RandomGenerator ranGen, std::vector<unsigned int> timeAtObservation,
		double startTime, int startDoy, double solarLatitude, double solarLongitude,
		double cellresolution){

	// create individual
	Individual *ind = new Individual(currentParameters, timeAtObservation,
			startTime, startDoy, solarLatitude, solarLongitude, cellresolution);
	// run steps for individual
	ind->move(iterations, maxDuration, env, ranGen);


	return ind;
	delete ind;
}

