/*
 * Solar.cpp
 *
 *  Created on: 25 Jan 2017
 *      Author: SH
 */

#include "Solar.h"

Solar::Solar(double latitude, double longitude) {

	sunrise = -999;
	sunset = -999;

	lat = latitude;		// lat = latitude in decimal degrees
	lon = longitude;	// lon = longitude in decimal degrees (negative == West)
	// TODO Auto-generated constructor stub
}

Solar::~Solar() {
	// TODO Auto-generated destructor stub
}

double Solar::deg2rad(double x){
	return(M_PI*x/180);
}

int Solar::sign(double x){
	if(x > 0) return 1;
	if(x == 0) return 0;
	if(x < 0) return -1;

	return(-999);
}

void Solar::solarcalc(int doy, double shiftSunrise, double shiftSunset){
	// doy = day of the year

	// This method is copied from:
	// https://www.r-bloggers.com/approximate-sunrise-and-sunset-times/
	// and they copied it from:
	// Teets, D.A. 2003. Predicting sunrise and sunset times.
	// The College Mathematics Journal 34(4):317-321.
	//
	// At the default location the estimates of sunrise and sunset are within
	// seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
	// with a mean of 2.4 minutes error.

	int earthr = 6378; // earth radius in Km

	double epsilon = deg2rad(23.45); // Radians between the xy-plane and the ecliptic plane

	double latRad = deg2rad(lat); // Convert lat to radians

	// Calculate offset of sunrise based on longitude (min)
	// If Long is negative, then the mod represents degrees West of
	// a standard time meridian, so timing of sunrise and sunset should
	// be made later.
	double timezone = -4.0 * fmod(fabs(lon), 15) * sign(lon);

	double r = 149598000; // The earth's mean distance from the sun (km)

	double theta = 2 * M_PI / 365.25 * (double) (doy-80);

	double zs = r * sin(theta) * sin(epsilon);
	double rp = sqrt(pow(r,2) - pow(zs,2));

	double t0 = 1440 / (2*M_PI)* acos((earthr - zs * sin(latRad)) / (rp * cos(latRad)));

	double that = t0 + 5; // kludge adjustment for the radius of the sun

	double n = 720 - 10 * sin(4*M_PI* (double) (doy-80) / 365.25) + 8 * sin(2*M_PI* (double) doy / 365.25); // Adjust "noon" for the fact that the earth's orbit is not circular:

	sunrise = (n - that + timezone + shiftSunrise) * 1; // in min
	sunset = (n + that + timezone + shiftSunset) * 1; // in min
}
