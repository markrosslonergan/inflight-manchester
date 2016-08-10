#ifndef CHANNEL_H_
#define CHANNEL_H_

#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "fourmomentum.h" // defines class fourmomentum
#include "sterile_flux.h" // defines class initial_sterile

#define CHAN_ELECPOSI 0
#define CHAN_ELECPI 1
#define CHAN_MUONPI 2
#define CHAN_NUPI0 3
#define CHAN_GAMMA 4
class initial_sterile;

typedef struct OBSERVABLES { // this is a struct of relevant observables for two (visible) body final states
	double E_sum; 	
	double Th_sum; 
	double AngSep; 
	double E_sterile; 
	double Th_sterile; 
	double E_high; 
	double Th_high; 
	double E_low; 
	double Th_low;
        std::vector<double > P_high;
	std::vector<double > P_low;
	double Minvar;
	double FS_AngSep; //The foreshortened angular separation.

	int chan_identifier;

	double E_high_smear;
	double E_low_smear;
	double Th_high_smear;
	double Th_low_smear;
	double AngSep_smear;
	double FS_AngSep_smear;
	double E_sum_smear;
	double Th_sum_smear;
	double E_sterile_smear;
	double Th_sterile_smear;
	double Minvar_smear;


	} OBSERVABLES;

class twoIP_channel { //This is the mother class for all decay channels (into two Ionising Particles)

public:
	twoIP_channel(gsl_rng * g, std::vector<double> input);

	fourmomentum IP1;	//first outgoing particle 4 momentum.
	fourmomentum IP2;	//second outgoing particle 4 momentum.
	int chan_identifier;
	gsl_rng * r;
	std::vector<double> model_params;

	int observables(OBSERVABLES * output, gsl_rng *g);
	virtual int decayfunction(initial_sterile nuS);	
};

/* ###############################
   
   Below here we have a derived class for each channel

   ############################### */

/* ########################################################################

	This is the nu_s \to \nu e+ e- channel (off-shell Zprime).

   ######################################################################## */

class threebody : public twoIP_channel {

public:
	threebody(gsl_rng * g, std::vector<double> input);
	int decayfunction(initial_sterile nuS);

	struct PDF_CHOICE { 
		double Enu; 
		double cosThnu; 
		double Phinu; };

//	typedef double (*threebody_pdf_function)(double, double, double, double, void *);

private:
	int computeLabFrameVariables(double mS, double Es, double costhS, double phiS, double restFrameParams[3]);
	double pdf_function(double x, double y, double mS, double mZprime, void * pointer);
//	struct PDF_CHOICE choose_from_pdf(gsl_rng * r, double mS, double mZprime, threebody_pdf_function pdf);
	int rotor(double theta, double phi, double vector[3]);
	int drawRestFrameDist(gsl_rng * r, double mS, double mZprime, double output[3]);
}; 


/* ########################################################################

	This is the nu_s \to \nu Zprime \to \nu e+ e- channel (on-shell Zprime).

   ######################################################################## */

class Zprimeresonance : public twoIP_channel {

public: 
	Zprimeresonance(gsl_rng * g, std::vector<double> input);
	int decayfunction(initial_sterile nuS);

private:
	double fourvec_costheta(double FOURVEC[4]);
	double fourvec_cosphi(double FOURVEC[4]);
	double rot_boost(double costh, double phi, double gam, double FOURVEC[4]);	
}; 

/* ########################################################################

	This is for a generic two body with a two (massive)  particle final state.

	This might be confusing... but when initialized you pass it two masses
	(ma and mb, in order).  The OBSERVABLES struct which is populated is then the
same as used for the eplus and eminus channel (defined at the top) -- but with
	E_high, costheta_high refering to the first particle (mass ma) and E_low ,
	costheta_low referring to the second particle (mass mb). It seems silly to
	actually order these different particles.

	Of course, we could come up with a more general OBSERVABLES function to
	make this a bit more logical.

   ######################################################################## */

class twobody : public twoIP_channel {

public: 
	twobody(gsl_rng * g, std::vector<double> input);
	int decayfunction(initial_sterile nuS);

}; 



#endif
