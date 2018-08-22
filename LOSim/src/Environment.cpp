/*
 * Environment.cpp
 *
 *  Created on: 14 Mar 2016
 *      Author: SH
 */

#include "Environment.h"

Environment::Environment(std::vector<std::vector<double> >*habitatSuitability //unsigned int res,
		//double xmin, double xmax, double ymin, double ymax
		) {
	// TODO Auto-generated constructor stub
	habitat = habitatSuitability;
//	extent[0] = xmin; // xmin, xmax, ymin, ymax
//	extent[1] = xmax;
//	extent[2] = ymin;
//	extent[3] = ymax;
//	resolution = res; // cell size, only quadratic
//	dim[0] =  (extent[1] - extent[0]) / resolution; // rows
//	dim[1] = (extent[3] - extent[2]) / resolution; // columns
	// TODO calculate DIMS
}

Environment::~Environment() {
	// TODO Auto-generated destructor stub
}

//double Environment::getExtent() const{
//	return extent;
//}


//void Environment::setExtent(double[] newExtent){
//	extent = newExtent;
//}
