/*
 *_________ _        _______  _       _________ _______          _________
 *\__   __/( (    /|(  ____ \( \      \__   __/(  ____ \|\     /|\__   __/
 *   ) (   |  \  ( || (    \/| (         ) (   | (    \/| )   ( |   ) (   
 *   | |   |   \ | || (__    | |         | |   | |      | (___) |   | |   
 *   | |   | (\ \) ||  __)   | |         | |   | | ____ |  ___  |   | |   
 *   | |   | | \   || (      | |         | |   | | \_  )| (   ) |   | |   
 *___) (___| )  \  || )      | (____/\___) (___| (___) || )   ( |   | |   
 *\_______/|/    )_)|/       (_______/\_______/(_______)|/     \|   )_(   
 *
 *		Inflight, Event generator for sterile decays at SBL facilities
 *
 *		If you have any questions, queries or comments please contact the authors;
 *			 mark.ross-lonergan@durham.ac.uk 
 *			 or 
 *			 peter.ballett@durham.ac.uk
 *
 *		The authors make no guarrentee of the behaviour, stability or bug-free-ness of this code.
 *		Use is at own risk.
 *
 */

#ifndef DETECTOR_H_
#define DETECTOR_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "channel.h"

#define ACCEPTED 0
#define REJECTED 1

#define DET_NOCUTS 0
#define DET_SBND 1
#define DET_MUBOONE 2
#define DET_ICARUS 3


double smear_energy(double En,double ms, double Percen, gsl_rng * r);
double smear_angle(double En, double Percen,gsl_rng * r);



struct OBSERVABLES;

class detector {

public: 

	virtual int accept(OBSERVABLES * Obs);
	int flag;
};


class nocuts : public detector {

public:
	int accept(OBSERVABLES * Obs);

};


class muBooNE : public detector {

public:
	muBooNE();
	int accept(OBSERVABLES * Obs);

	double Energy_threshold; //Below this events aren't registered.

	double AngSep_threshold; //Below this electorn positron pairs are seen as a single track.

	double Energy_ratio_threshold; 	// If E_low/E_high is below this the 
					// event is treated as a single track.

};

class muBooNE_ee : public detector {

public:
	muBooNE_ee();
	int accept(OBSERVABLES * Obs);

	double AngSep_threshold_low; //Below this electorn positron pairs are seen as a single track.
	double AngSep_threshold_high; //background cut on events that are seperated by more than this
	double Energy_sum_threshold; //onyl accept energetic enough events

					// event is treated as a single track.

};

class muBooNE_gamma : public detector {

public:
	muBooNE_gamma();
	int accept(OBSERVABLES * Obs);

	double Energy_threshold; //Below this events aren't registered.
	double AngSep_threshold; //Below this electorn positron pairs are seen as a single track.
	double Ang_threshold; //only accept events that are going very very forward


};


class muBooNE_epi : public detector {

public:
	muBooNE_epi();
	int accept(OBSERVABLES * Obs);

	double Energy_sum_threshold; //Below this events aren't registered.
	double AngSep_threshold; //Below this electorn positron pairs are seen as a single track.


};

class muBooNE_mupi : public detector {

public:
	muBooNE_mupi();
	int accept(OBSERVABLES * Obs);
	double Energy_sum_threshold; //Below this events aren't registered.
	double AngSep_threshold; //Below this electorn positron pairs are seen as a single track.


};

class muBooNE_nupi0 : public detector {

public:
	muBooNE_nupi0();
	int accept(OBSERVABLES * Obs);
	double Energy_reco_threshold; //Below this events aren't registered.
	double Ang_reco_threshold; //Below this electorn positron pairs are seen as a single track.


};


#endif
