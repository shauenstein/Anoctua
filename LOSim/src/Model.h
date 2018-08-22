/*
 * Model.h
 *
 *  Created on: 14 Mar 2016
 *      Author: SH
 */

#ifndef MODEL_H_INCLUDED__
#define MODEL_H_INCLUDED__
#include "Environment.h"
#include "Individual.h"
#include "RandomGenerator.h"

//#include <iostream>

class Model {
public:
	Model(std::vector<std::vector<double> >*habitatSuitability);
	virtual ~Model();
	Individual* runSimulation(std::vector<double> currentParameters,
			unsigned int iterations, int maximumDuration,
			RandomGenerator ranGen, std::vector<unsigned int> timeAtObservation,
			double startTime, int startDoy, double solarLatitude, double solarLongitude,
			double cellresolution);
	Environment* env;
};

#endif /* MODEL_H_ */
