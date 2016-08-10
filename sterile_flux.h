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
	std::vector<double > fluxlist;
	std::vector<double > elist;

	double fmax;
	double mass;

	double get_event(gsl_rng *r);


};


#endif
