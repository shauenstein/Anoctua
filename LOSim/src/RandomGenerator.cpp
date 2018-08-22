/*
 * RandomGenerator.cpp
 *
 *  Created on: 28 Apr 2016
 *      Author: SH
 */

#include "RandomGenerator.h"

RandomGenerator::RandomGenerator() {
	// TODO Auto-generated constructor stub

}

RandomGenerator::~RandomGenerator() {
	// TODO Auto-generated destructor stub
}


std::default_random_engine & RandomGenerator::my_engine( )
{
 static std::default_random_engine e;
 return e;
}

void RandomGenerator::seedrand(unsigned int s) { my_engine().seed(s); }

double RandomGenerator::randomdouble(double min, double max)
{
	 std::uniform_real_distribution<double> unif(min,max);
	 return unif(my_engine());
}

double RandomGenerator::rGaussian(double mean, double sd)
{
	std::normal_distribution<double> norm(mean, sd);

	return norm(my_engine());
}

int RandomGenerator::rPoisson(double lambda)
{
	std::poisson_distribution<int> pois(lambda);

	return pois(my_engine());
}

double RandomGenerator::rLogNormal(double logMean, double logSd)
{
	std::lognormal_distribution<double> lnorm(logMean, logSd);

	return lnorm(my_engine());
}

double RandomGenerator::rLogNormalTruncated(double logMean, double logSd, double maximum)
{
	double value = maximum + 1;
	while(value > maximum){
		value = rLogNormal(logMean, logSd);
	}
	return value;
}

int RandomGenerator::rBinomial(double prob)
{
	std::binomial_distribution<int> rbinom(1, prob);

	return rbinom(my_engine());
}

int RandomGenerator::multinomialDraw(double *chances, int size, double max)
{

	double min = 0.0;
	double x = randomdouble(min, max);
	//std::cout << x << std::endl;

	bool go = true;
	int start = 0;
	int end = size;
	int mid = size/2;
	double midValue = chances[mid];
	while(go)
	{
		if(x < midValue)
		{
			end = mid;
			mid = (start+end)/2;
			midValue = chances[mid];
		}else{
			start = mid;
			mid = (start+end)/2;
			midValue = chances[mid];
		}

		if (start == mid){
			if (mid == 0 && x < midValue ) return 0;
			else return mid + 1;
		}
	}
	return -999;
}

