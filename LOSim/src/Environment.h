/*
 * Environment.h
 *
 *  Created on: 14 Mar 2016
 *      Author: SH
 */

#ifndef ENVIRONMENT_H_INCLUDED__
#define ENVIRONMENT_H_INCLUDED__
#include <vector>
//#include <boost/numeric/ublas/storage.hpp>

class Environment{
public:
	// default constructor
	//Environment();

	// overload constructor
	Environment(std::vector<std::vector<double> >*habitatSuitability //unsigned int res,
			//double xmin, double xmax, double ymin, double ymax
			);

	// destructor
	virtual ~Environment();

//private:
	// accesor functions // just as an example, only needed for private variables
	//double getExtent() const;

	// mutator functions
	//void setExtent(double[])

	std::vector<std::vector<double> >* habitat;
//	double extent[4];
//	unsigned int resolution;
//	unsigned int dim[2];
};
#endif /* ENVIRONMENT_H_ */
