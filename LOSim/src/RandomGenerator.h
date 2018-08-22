/*
 * RandomGenerator.h
 *
 *  Created on: 28 Apr 2016
 *      Author: SH
 */

#ifndef RANDOMGENERATOR_H_INCLUDED__
#define RANDOMGENERATOR_H_INCLUDED__

#include <random>
#include <vector>
//#include <iostream>

class RandomGenerator {
public:
	RandomGenerator();
	virtual ~RandomGenerator();

	int multinomialDraw(double *chances, int size, double max);
	std::default_random_engine & my_engine();
	void seedrand(unsigned int s);
	double randomdouble(double min, double max);
	double rGaussian(double mean, double sd);
	int rPoisson(double lambda);
	int rBinomial(double prob);
	double rLogNormal(double logMean, double logSd);
	double rLogNormalTruncated(double logMean, double logSd, double maximum);
};

#endif /* RANDOMGENERATOR_H_ */
