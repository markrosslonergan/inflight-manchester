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

#ifndef STERILE_FLUX_H_
#define STERILE_FLUX_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>


#include "fourmomentum.h"
#include "detector.h"

class initial_sterile {

public:
	double mass;
	double energy;
	double costhS;
	double phiS;
	fourmomentum labframeP;

	initial_sterile(double M, double E, double in_costhS, double in_phiS);
};

double getEvents(double mS, int detector, double events[][2]);


double interpolate(std::vector<double> xlist, std::vector<double> weightlist, double xin);

struct fluxfile
{
	fluxfile(std::string name, double min);
	double getFlux(double Ein);

	std::string filename;
	std::vector<double> fluxlist;
	std::vector<double> elist;

	double fmax;
	double mass;

	double get_event(gsl_rng *r);


};


#endif
