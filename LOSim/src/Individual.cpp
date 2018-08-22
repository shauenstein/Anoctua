/*
 * Individual.cpp
 *
 *  Created on: 31 Mar 2016
 *      Author: SH
 */

#include "Individual.h"

Individual::Individual(std::vector<double> currentParameters,
		std::vector<unsigned int> Observationtime,
		double startT,
		int startD,
		double solarLat,
		double solarLon,
		double res) {

	cellres = res;
	sol = new Solar(solarLat, solarLon);
	tObs = Observationtime;

	// start coordinates
	x.push_back((int) currentParameters.at(0));
	y.push_back((int) currentParameters.at(1));
	turningAngle.push_back(NAN);

	stepDistance.push_back(0.0);
	energy = 1;
	energystamp.push_back(energy);

	// day and time
	daytime  = startT * 60; // to min
	doy = startD;

	active = true; // activity status // initialise

	// Parameters
	habitatPreference =	exp(currentParameters.at(2))-1;

	stepShape = currentParameters.at(3);
	stepScale = currentParameters.at(4);

	directionalBias = 1 / currentParameters.at(5);

	DispersalRestingTimeLogMean = currentParameters.at(6); // min
	DispersalRestingTimeLogSd = currentParameters.at(7); // min

	maximumEffort = currentParameters.at(8);

	roostLambda = currentParameters.at(9);

	if(!tObs.empty()){
		observationModel = true;
		observationError = currentParameters.at(10);
	    rscRange = (int) currentParameters.at(11);
		numberOfSteps = 0;
		startDoy = startD; // copy for observation model
		startDaytime = startT * 60; // copy for observation model
	}else{
		observationModel = false;
	}

	validityCode = 0;


	// TODO Auto-generated constructor stub

}

Individual::~Individual() {
	// TODO Auto-generated destructor stub
}

// move (steps)
void Individual::move(unsigned int nSteps, int maximumDuration, Environment *env, RandomGenerator ranGen){

	// initiate gamma step length distribution
	boost::math::gamma_distribution<double> stepDistr(stepShape, stepScale);
	int perceptionRange = (int) ceil(quantile(stepDistr, 0.95));


	// initiate gaussian turning angle distribution
	boost::math::normal_distribution<double> angleDistr(0, directionalBias);

	// set maximumDuration if observationModel
	if(observationModel){
		maximumDuration = tObs.at(tObs.size()-1);
	}

	// scanning:
	// create array with Number of vertical cells along horizontal axis

	int NumberOfCells = 0;
	int oppositeLegs[2 * perceptionRange + 1]; // array of opposite legs
	int subscript = 0; // loop subscript
	for(int pR = - perceptionRange; pR <= perceptionRange; pR++, subscript++){
		// get indices for area within perception
		oppositeLegs[subscript] = floor(sqrt(pow(perceptionRange, 2.0) - pow(pR, 2.0)));
		NumberOfCells += 2 * oppositeLegs[subscript] + 1;
	}

	// fill in starting conditions

	// calculate sunrise and sunset times
	sol->solarcalc(doy, -25.0, 25.0);

	// evaluate starting activity status
	updateActivity();

	// additional start point states
	//envVal.push_back(mooreEnv(env, x.at(0), y.at(0)));
	envVal.push_back(env->habitat->at(y.at(0)).at(x.at(0)));

	if(!active){
		// only day rest
		activitystatus.push_back(1);
		active = true;
		timestamp.push_back(sol->sunset - daytime);
		daytime = sol->sunset;
	}else{
		activitystatus.push_back(0);
		timestamp.push_back(0.0);
	}
	time.push_back(daytime);

	// start moving
	//
	for(unsigned int iter = 0; iter < nSteps; iter++){
		//

		// check if still within range of environment
		if((x.at(iter) - perceptionRange) <= 2 ||
			(x.at(iter) + perceptionRange) >= (int) (env->habitat->at(0).size() - 1) ||
			(y.at(iter) - perceptionRange) <= 2 ||
			(y.at(iter) + perceptionRange) >= (int) (env->habitat->size() - 1)){
			//
			validityCode = 1;
			numberOfSteps = iter;
			// break out of iteration loop
			iter = nSteps;
		}else{
			// compute weights
			double cumulativeWeights[NumberOfCells];
			//
			double sumOfWeights = 0.0;
			//
			double cumEnv = 0.0;
			//
			subscript = 0; // reset loop subscript
			int currentCell = 0; // Cell subscript in cumulativeWeights array
			//
			// horizontal loop, x-direction
			for(int h = - perceptionRange; h <= perceptionRange; h++, subscript++){
				// vertical loop, y-direction
				for(int v = oppositeLegs[subscript]; v >= -oppositeLegs[subscript]; v--, currentCell++){

					// sum up habitat suitability
					cumEnv += env->habitat->at(y.at(iter) + v).at(x.at(iter) + h);

					// skip point of origin, give a weight of zero
					if(h == 0 && v == 0){
						sumOfWeights += 0;

					}else{
						// step after being inactive (or first step)
						if(activitystatus.at(iter) != 0 || iter == 0){
							sumOfWeights +=
								//pow(mooreEnv(env, x.at(iter) + v, y.at(iter) + h), habitatPreference) * // habitat prefence
								pow(env->habitat->at(y.at(iter) + v).at(x.at(iter) + h), habitatPreference) *
								pdf(stepDistr, sqrt(pow(h, 2) + pow(v, 2))); // step length pdf from gamma

						}else{ // step while being active (dispersing)
							// current angle
							double preVec[2] = {(double) (x.at(iter) - x.at(iter - 1 )), (double) (y.at(iter) - y.at(iter - 1))};
							double currentAngle = ((double) h * preVec[0] + (double) v * preVec[1]) /
									(sqrt(pow((double) h, 2) + pow((double) v, 2)) * sqrt(pow(preVec[0], 2) + pow(preVec[1], 2)));
							currentAngle = acos(currentAngle);
							// problem when directly turning back i.e. acos(-1) makes NAN in c##
							if(std::isnan(currentAngle)){
								currentAngle = M_PI;
							}

							sumOfWeights += //pow(mooreEnv(env, x.at(iter) + v, y.at(iter) + h), habitatPreference) * // habitat prefence
									pow(env->habitat->at(y.at(iter) + v).at(x.at(iter) + h), habitatPreference) *
								pdf(stepDistr, sqrt(pow(h, 2) + pow(v, 2))) * // step length pdf from gamma
								pdf(angleDistr, currentAngle);
						}
					}
					cumulativeWeights[currentCell] = sumOfWeights;
				}
			}

			// standardise cumulated habitat suitability
			cumEnv = cumEnv / NumberOfCells;

			// draw from multinomial distribution
			int selected = ranGen.multinomialDraw(cumulativeWeights, NumberOfCells, sumOfWeights);
			while(selected == -999){
				selected = ranGen.multinomialDraw(cumulativeWeights, NumberOfCells, sumOfWeights);
			}

			// find selected x coordinate
			int relativeX = 0;

			subscript = 0; // reset loop subscript
			for(int pR = - perceptionRange; pR <= perceptionRange; pR++, subscript++){
				selected -= 2 * oppositeLegs[subscript] + 1;
				if(selected < 0){
					relativeX = pR;
					break;
				}
			}

			// find selected y coordinate
			int relativeY = 0;
			int LegArray[2 * oppositeLegs[subscript] + 1];

			int subscript2 = 0; // reset loop subscript
			for(int i = oppositeLegs[subscript]; i >= -oppositeLegs[subscript]; i--, subscript2++){
				LegArray[subscript2] = i;
			}

			relativeY = LegArray[(2 * oppositeLegs[subscript]) + 1 + selected];

			// fill in conditions at iter + 1

			// write into coordinates matrix
			x.push_back(x.at(iter) + relativeX);
			y.push_back(y.at(iter) + relativeY);
			// write environmental value into coordinates matrix
			//envVal.push_back(mooreEnv(env, x.at(iter + 1), y.at(iter + 1)));
			envVal.push_back(env->habitat->at(y.at(iter + 1)).at(x.at(iter + 1)));

			// update step distance
			stepDistance.push_back(sqrt(pow(relativeX, 2) +
										pow(relativeY, 2)));

			// reduce energy state, if active
			if(active){
				energy -= (stepDistance.at(iter + 1) / maximumEffort);
				checkEnergyLimit();
			}

			// turning Angle
			if((iter) > 0){ // turning angle from third step on in hindsight
				double preVec[2] = {(double) (x.at(iter) - x.at(iter - 1)), (double) (y.at(iter) - y.at(iter - 1))};
				double curVec[2] = {(double) (x.at(iter + 1) - x.at(iter)), (double) (y.at(iter + 1) - y.at(iter))};

				turningAngle.push_back(acos((curVec[0] * preVec[0] + curVec[1] * preVec[1]) /
											(sqrt(pow(curVec[0], 2) + pow(curVec[1], 2)) *
													sqrt(pow(preVec[0], 2) + pow(preVec[1], 2)))));
			}else{ // for first step no angle calculation possible
				turningAngle.push_back(NAN);
			}

			// check activitystatus
			if(!active){
				// evaluate if day rest or real roosting
				if(ranGen.rBinomial((cumEnv + 1-energy)/2)){ // real roosting
					activitystatus.push_back(2 + ranGen.rPoisson(roostLambda));
					doy += activitystatus.at(iter + 1) - 1;
					sol->solarcalc(doy, -25.0, 25.0);
					timestamp.push_back(timestamp.at(iter) +
							(activitystatus.at(iter + 1) - 1)*24*60 +
							sol->sunset - daytime +
							ranGen.rLogNormalTruncated(DispersalRestingTimeLogMean, DispersalRestingTimeLogSd, 30));
				}else{ // day rest
					activitystatus.push_back(1);
					timestamp.push_back(timestamp.at(iter) + sol->sunset - daytime +
							ranGen.rLogNormalTruncated(DispersalRestingTimeLogMean, DispersalRestingTimeLogSd, 30));
				}
				// linear energy uptake as a function of the habitat quality and the time
				// full energy after 3 days when energy = 0 and habitat = 1
				energy += cumEnv * (timestamp.at(iter + 1) - timestamp.at(iter)) / (3*24*60);
				checkEnergyLimit();
				active = true;
				daytime = sol->sunset;
			}else{
				timestamp.push_back(timestamp.at(iter) +
						ranGen.rLogNormalTruncated(DispersalRestingTimeLogMean, DispersalRestingTimeLogSd, 30));
				daytime += timestamp.at(iter + 1) - timestamp.at(iter);
				if(daytime > 24*60){
					doy += 1;
					daytime -= 24*60;
					sol->solarcalc(doy, -25.0, 25.0);
				}
				// evaluate and update activity status
				activitystatus.push_back(0);
				updateActivity();
			}

			// update energystamp and time
			energystamp.push_back(energy); // energy
			time.push_back(daytime);

			// check validity
			if(iter+1 == nSteps){ // iteration loop ends, irrespective of observation model
				if(observationModel){ // observations not fully consumed
					validityCode = 2;
				}
				numberOfSteps = iter + 1;
			}
			if(maximumDuration < timestamp.at(iter +1)){
				numberOfSteps = iter + 1;
				iter = nSteps; // end loop
			}
		// end range check
		}
	// end iteration loop
	}
	if(observationModel){
		// run Observation model
		observe(ranGen, env);
	}
// end final
}

void Individual::observe(RandomGenerator ranGen, Environment *env){
	// states at position 0 for both simulated and observed time
	xObs.push_back(x.at(0));
	yObs.push_back(y.at(0));
	envValObs.push_back(envVal.at(0));
	stepDistanceObs.push_back(0.0);
	turningAngleObs.push_back(NAN);
	rsc.push_back(NAN);

	unsigned int match = 1; // loop: start check at position 1 in simulated time
	int currentDoy;
	double currentDaytime;

	// states along tObs
	for(unsigned int obs = 1; obs < tObs.size(); obs++){
		// match observed and simulated time
		match = locate(tObs.at(obs), match);
		// get states at match
		if(activitystatus.at(match) < 2){ // day rest or active dispersal
			xObs.push_back(x.at(match) + ranGen.rGaussian(0, observationError));
			yObs.push_back(y.at(match) + ranGen.rGaussian(0, observationError));
		}else{ // if real roosting
			// calculate daytime and doy
			currentDoy = startDoy + (int)(tObs.at(obs)/(24*60));
			currentDaytime = startDaytime + (tObs.at(obs) - (currentDoy - startDoy)*24*60);
			if(currentDaytime > 24*60){
				currentDaytime -= 24*60;
				currentDoy += 1;
			}
			// update solar information
			sol->solarcalc(currentDoy, -25.0, 25.0);

			if(currentDaytime > sol->sunrise && currentDaytime < sol->sunset){ // day rest
					xObs.push_back(x.at(match) + ranGen.rGaussian(0, observationError));
					yObs.push_back(y.at(match) + ranGen.rGaussian(0, observationError));
				}else{ // night roost
					xObs.push_back(x.at(match) + ranGen.rGaussian(0, observationError + 50/cellres));
					yObs.push_back(y.at(match) + ranGen.rGaussian(0, observationError + 50/cellres));
				}
		}
		// observed environment
		if(xObs.at(obs)  <= 2 || // range check
			 xObs.at(obs) >= (int) (env->habitat->at(0).size() - 1) ||
			 yObs.at(obs) <= 2 ||
			 yObs.at(obs) >= (int) (env->habitat->size() - 1)){
			envValObs.push_back(NAN);
		}else{
			//envValObs.push_back(mooreEnv(env, xObs.at(obs), yObs.at(obs)));
			envValObs.push_back(env->habitat->at(yObs.at(obs)).at(xObs.at(obs)));
		}
		// resource selection coefficient
		if(xObs.at(obs)  <= rscRange+1 || // range check
		   xObs.at(obs) >= (int) (env->habitat->at(0).size() - rscRange) ||
		   yObs.at(obs) <= rscRange+1 ||
		   yObs.at(obs) >= (int) (env->habitat->size() - rscRange)){
		  rsc.push_back(NAN);
		}else{
		  rsc.push_back(rsc_calc(env, xObs.at(obs), yObs.at(obs)));
		}
		// step distance
		stepDistanceObs.push_back(sqrt(pow(xObs.at(obs) - xObs.at(obs-1), 2) +
										pow(yObs.at(obs) - yObs.at(obs-1), 2)));
		// turning Angle
		if((obs) > 1){ // turning angle from third step on in hindsight
			double preVec[2] = {(double) (xObs.at(obs-1) - xObs.at(obs - 2)), (double) (yObs.at(obs-1) - yObs.at(obs - 2))};
			double curVec[2] = {(double) (xObs.at(obs) - xObs.at(obs-1)), (double) (yObs.at(obs) - yObs.at(obs-1))};
			double TAO = acos((curVec[0] * preVec[0] + curVec[1] * preVec[1]) /
					(sqrt(pow(curVec[0], 2) + pow(curVec[1], 2)) *
					sqrt(pow(preVec[0], 2) + pow(preVec[1], 2))));
			// problem when directly turning back i.e. acos(-1) makes NAN in c##
			if(std::isnan(TAO)){
				TAO = M_PI;
			}
			turningAngleObs.push_back(TAO);
		}else{ // for first step no angle calculation possible
			turningAngleObs.push_back(NAN);
		}

		// for validity code = 2, not entirely consumed observations
		if(match == numberOfSteps - 1){
			obs = tObs.size(); // end loop
		}
	}
}

// locate position where observed time matches simulated time
unsigned int Individual::locate(unsigned int targetTime, unsigned int startPosition){
	int position;
	for(unsigned int l = startPosition; l <= numberOfSteps; l++){
		if(timestamp.at(l) >= targetTime){
			position = l;
			l = numberOfSteps; // end loop
		}
	}
	return position;
}

// moore neighbourhood
double Individual::mooreEnv(Environment *envir, int xcoord, int ycoord){
	return (envir->habitat->at(ycoord).at(xcoord) + // mid
			envir->habitat->at(ycoord + 1).at(xcoord - 1) + // bottom left
			envir->habitat->at(ycoord + 1).at(xcoord) + // bottom
			envir->habitat->at(ycoord + 1).at(xcoord + 1) + //bottom right
			envir->habitat->at(ycoord).at(xcoord + 1) + // right
			envir->habitat->at(ycoord - 1).at(xcoord + 1) + // top right
			envir->habitat->at(ycoord - 1).at(xcoord) + // top
			envir->habitat->at(ycoord - 1).at(xcoord - 1) + // top left
			envir->habitat->at(ycoord).at(xcoord - 1)) // left
			/ 9;
}

double Individual::rsc_calc(Environment *envir, int xcoord, int ycoord){
  double out;
//  double minEnv = 1;
//  double maxEnv = 0;
  int sEnv = 1; // number of env. values smaller than at current location
  int lEnv = 1; // number of env. values larger than at current location
  double curEnv = envir->habitat->at(ycoord).at(xcoord); // current env. value

  int rscLegs[2 *  rscRange + 1]; // array of opposite legs
  int rscIndex = 0; // loop subscript
  for(int aR = - rscRange; aR <= rscRange; aR++, rscIndex++){
    // get indices for area within perception
    rscLegs[rscIndex] = floor(sqrt(pow(rscRange, 2.0) - pow(aR, 2.0)));
  }

  // reset subscript
  rscIndex = 0;

  // loop through cells within rscRange
  for(int h = - rscRange; h <= rscRange; h++, rscIndex++){
    // vertical loop, y-direction
    for(int v = rscLegs[rscIndex]; v >= -rscLegs[rscIndex]; v--){
      // skip point of origin, give a weight of zero
      if(h != 0 && v != 0){
        // check if env is larger or smaller than curEnv
        if(envir->habitat->at(ycoord + v).at(xcoord + h) > curEnv){
          lEnv += 1;
        }else{
          if(envir->habitat->at(ycoord + v).at(xcoord + h) < curEnv){
            sEnv += 1;
          }
        }
//        // update min and max env. values
//        if(envir->habitat->at(ycoord + v).at(xcoord + h) < minEnv){
//          minEnv = envir->habitat->at(ycoord + v).at(xcoord + h);
//        }
//        if(envir->habitat->at(ycoord + v).at(xcoord + h) > maxEnv){
//          maxEnv = envir->habitat->at(ycoord + v).at(xcoord + h);
//        }
      }
    }
  }
  // ratio divided by range
//  if((maxEnv - minEnv) == 0){
//    out = NAN;
//  }else{
//    out = (lEnv/ (double) sEnv) / (maxEnv - minEnv);
  out = lEnv/ (double) sEnv;
//  }

  return out;
}

// update activity
void Individual::updateActivity(){
	// evaluate starting activity status
	if(daytime > sol->sunrise && daytime < sol->sunset){
		active = false;
	}else{
		active = true;
	}
}

void Individual::checkEnergyLimit(){
	if(energy < 0) energy = 0; // cannot be less than zero
	if(energy > 1) energy = 1; // cannot be more than one
}
