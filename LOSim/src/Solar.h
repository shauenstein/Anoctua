/*
 * Solar.h
 *
 *  Created on: 25 Jan 2017
 *      Author: SH
 */

#ifndef SOLAR_H_INCLUDED__
#define SOLAR_H_INCLUDED__

#include <math.h>
#include <vector>
#include "RandomGenerator.h"

// define pi which was dropped in C99
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

class Solar {
public:
	Solar(double latitude, double longitude);
	virtual ~Solar();

	double sunrise;
	double sunset;
	double lat;
	double lon;

	void solarcalc(int doy, double shiftSunrise, double shiftSunset);
	double deg2rad(double x);
	int sign(double x);

};

#endif /* SOLAR_H_ */
