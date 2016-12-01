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
 * Copyright (2016) 
 * Mark Ross-Lonergan mark.ross-lonergan@durham.ac.uk 
 * Peter Ballett peter.ballett@durham.ac.uk
 */
#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define no_argument 0
#define required_argument 1
#define optional_argument 2

#include "fourmomentum.h" // Defines a class for fourmomenta; makes getting directions and 3-momenta (and
			  // all that) easier.

#include "sterile_flux.h" // Defines getEvents, and describes a class whose objects denote a single 
			  // incoming sterile (mass, fourmomentum).

#include "channel.h"	  // Includes the mother class for two body final state
			  // decays, and the derived classes for the two channels we've studied so far.

#include "detector.h"	  // This includes the detector-specific cut functions.


#define MPION  0.13957
#define MPI0   0.13497
#define MKAON  0.49367
#define MMU    0.10566
#define	ME     0.00051



int output_distributions(gsl_rng * r, detector * DETECTOR, twoIP_channel * CHAN, double mS, bool hepornot, std::string flux_name, int NUMEVENTS);
int migration_matrix(detector * DETECTOR, twoIP_channel * CHAN, double mS, gsl_rng *g);	
int efficiency_matrix(detector * DETECTOR, twoIP_channel * CHAN, double mS, gsl_rng *g);	

/* ########## Main function ############### */
int main(int argc, char* argv[])
{

const gsl_rng_type * T; // Standard invocation prayer to the RNG
gsl_rng * r;
gsl_rng_env_setup();

T = gsl_rng_default;
r = gsl_rng_alloc (T);

double mS = 0.1500; 	 		// These are the default values if no command line parameters are defined. 
int channel_flag = CHAN_ELECPI; 	// This one too.
int detector_flag = DET_NOCUTS;	 	// This one too.
int matrix_flag = 0; 	 		// This one too.
int eff_flag = 0; 			// This one too.
bool hepevt = false;
std::string flux_name = "flux.dat";

int NUMEVENTS = 100;

int c;
opterr = 0;

int index; 
int iarg = 0;
opterr=1;
const struct option longopts[] = 
{
	{"mass",	 	required_argument, 	0, 'm'},
	{"migration",	 	no_argument, 		0, 'M'},
	{"efficiency", 		no_argument, 		0, 'E'},
	{"flux-file",		required_argument,	0, 'f'},
	{"help",		no_argument,		0, 'h'},
	{"hepevt", 		no_argument, 		0, 'v'},
	{"using-no-detector",	no_argument, 		0, 'o'},
	{"using-muboone",	no_argument, 		0, 'b'},
	{"number",		required_argument,	0, 'n'},
	{"channel",		required_argument,	0, 'C'},
	{0,			no_argument, 		0,  0},
};


while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "m:MEhvobC:n:f:", longopts, &index);

	switch(iarg)
	{
		case 'C':
			channel_flag = strtod(optarg,NULL);
			break;
		case 'm':
			mS = strtof(optarg,NULL);
			break;
		case 'v':
			hepevt= true;
			break;	
		case 'f':
			flux_name =optarg;
			break;
		//case 'M':
		//	matrix_flag = 1;
		//	break;
		//case 'E':
		//	matrix_flag = 1;
		//	eff_flag = 1;		
		//	break;
		case 'n':
			NUMEVENTS = strtod(optarg,NULL);
			break;
		case 'o':
			detector_flag = DET_NOCUTS;
			break;
		case 'b':
			detector_flag = DET_MUBOONE;
			break;
		case '?':
		std::cout<<"Abandon hope all ye who enter this value. "<<std::endl;
		case 'h':
		std::cout<<"\nInflight, an event generator for heavy sterile decays at SBL facilities"<<std::endl;
		std::cout<<"If you have any questions please contact the authors mark.ross-lonergan@durham.ac.uk or peter.ballett@durham.ac.uk"<<std::endl;
		std::cout<<"The authors make no guarrentee of this code"<<std::endl;
		std::cout<<"******************************************"<<std::endl;
		std::cout<<"Allowed arguments:"<<std::endl;
		std::cout<<"\t-m, --mass\t\t\tSets the parent sterile mass. [default = 0.1500]"<<std::endl;
		std::cout<<"\t-n, --number\t\t\tHow many events will we generate? [Default 100]"<<std::endl;
		std::cout<<"\t-C, --channel\t\t\tsets the decay channel [default = 1] \n\t\t\t\t\t\t0: 3-body (nu e e).\n\t\t\t\t\t\t1: Isotroic 2-body e pi\n\t\t\t\t\t\t2: Isotropic 2-body mu pi\n\t\t\t\t\t\t3: Isotropic 2-body nu Pi0\n\t\t\t\t\t\t4: Isotropic 2-body nu gamma\n\t\t\t\t\t\t5: 3-body nu mu mu\n\t\t\t\t\t\t6: 3-body nu mu e"<<std::endl;
		//std::cout<<"\t-M or --migration\t\t\tproduces a migration matrix for oberservable 'n' (argument probably doesn't work as yet)"<<std::endl;
		//std::cout<<"\t-E of --efficiency\t\t\tproduces an efficiency matrix for total energy. "<<std::endl;
		std::cout<<"\t--using-no-detector\t\truns with no detector cuts. [DEFAULT] "<<std::endl;
		std::cout<<"\t--using-muboone\t\t\truns with muBooNE detector cuts."<<std::endl;
		std::cout<<"\t-f, --flux-file\t\t\tThe name of the file containing flux to sample from, no normalisation needed. 1st column energy in GeV, 2nd column flux."<<std::endl;
		std::cout<<"\t-v, --hepevt\t\t\tOutputs to hepevt format."<<std::endl;
		std::cout<<"\t-h, --help\t\t\tDisplays this message"<<std::endl;
		std::cout<<"******************************************"<<std::endl;
		std::cout<<"Example usage, if file flux_150.dat exists in the directory then to generate 2000 N->mu-pi+ decays in the hepevt format either"<<std::endl;
		std::cout<<"  ./inflight -m 0.150 -n 2000 -C 2 -v -f flux_150.dat "<<std::endl;
		std::cout<<"  ./inflight --mass 0.150 --number 2000 --channel 2 --hepevt --flux-file flux_150.dat\n\n"<<std::endl;

			return 0;
	}

}




//Now we set up the decay channel object.
std::vector<double> model_params; //This should include any theoretical parameters which the model needs to know.

twoIP_channel *CHAN;

switch(channel_flag)
{
	case CHAN_ELECPOSI:
		model_params.push_back(91.19); //mediator mass
		model_params.push_back((double) CHAN_ELECPOSI); 
		model_params.push_back((double) CHAN_ELECPOSI); 
		CHAN = new threebody(r,model_params); 
//		std::cout<<"channel: e e"<<std::endl;
		break;
	case CHAN_ELECPI:
		model_params.push_back(ME); // the electron with an incorrect mass (to try and avoid spacelike vectors...).
		model_params.push_back(MPION); // the pion. 0.135
		model_params.push_back((double) CHAN_ELECPI); 
		CHAN = new twobody(r,model_params); 
//		std::cout<<"channel: e pi"<<std::endl;
		break;

	case CHAN_MUONPI:
		model_params.push_back(MMU); // the muon.
		model_params.push_back(MPION); // the pion.
		model_params.push_back((double) CHAN_MUONPI); 
		CHAN = new twobody(r,model_params); 
//		std::cout<<"channel: mu pi"<<std::endl;
		break;
	case CHAN_NUPI0:
		model_params.push_back(0.00); // the neutrino
		model_params.push_back(MPI0); // the pion0
		model_params.push_back((double) CHAN_NUPI0);
		CHAN = new twobody(r,model_params); 
//		std::cout<<"channel: mu pi"<<std::endl;
		break;
	case CHAN_GAMMA:
		model_params.push_back(91.19); //mediator mass
		model_params.push_back((double) CHAN_GAMMA); // the pion. 
		model_params.push_back((double) CHAN_GAMMA); // the pion. 
		CHAN = new threebody(r,model_params); 
//		std::cout<<"channel: e e"<<std::endl;
		break;
	case CHAN_MUMU:
		model_params.push_back(91.19); //mediator mass
		model_params.push_back((double) CHAN_MUMU); // the pion. 
		model_params.push_back((double) CHAN_MUMU); // the pion. 
		CHAN = new threebody(r,model_params); 
		break;
	case CHAN_MUE:
		model_params.push_back(91.19); //mediator mass
		model_params.push_back((double) CHAN_MUE); // the pion. 
		model_params.push_back((double) CHAN_MUE); // the pion. 
		CHAN = new threebody(r,model_params); 
		break;



	default:
		std::cout<<"ERROR: Bad channel specifier."<<std::endl;
		return -1;
}

//We define the detector cuts we would like to use. 
//
detector * DETECTOR;

switch(detector_flag)
{
	case DET_NOCUTS:
		DETECTOR = new nocuts(); 	// pseudo-detector that just allows every event.
//		std::cout<<"detector: no cuts"<<std::endl;
		break;
	case DET_MUBOONE:
		switch(channel_flag)
		{
		case CHAN_ELECPOSI:	
			DETECTOR = new muBooNE_ee(); 	// electron positron by off-shell Z-like.
//			std::cout<<"detector: muB ee"<<std::endl;
			break; 
		case CHAN_ELECPI:
			DETECTOR = new muBooNE_epi(); 	// isotropic electron pion.
//			std::cout<<"detector: muB e pi"<<std::endl;
			break; 
		case CHAN_MUONPI:
			DETECTOR = new muBooNE_mupi(); 	// isotropic muon pion.
//			std::cout<<"detector: muB mu pi"<<std::endl;
			break;
		case CHAN_NUPI0:
			DETECTOR = new muBooNE_nupi0(); 	// isotropic muon pion.
//			std::cout<<"detector: muB mu pi"<<std::endl;
			break;
		case CHAN_GAMMA:
			DETECTOR = new muBooNE_gamma(); 	// isotropic muon pion.
//			std::cout<<"detector: muB mu pi"<<std::endl;
			break;
		case CHAN_MUMU:
			DETECTOR = new muBooNE_ee(); 	// isotropic muon pion.
//			std::cout<<"detector: muB mu pi"<<std::endl;
			break;
		case CHAN_MUE:
			DETECTOR = new muBooNE_ee(); 	// isotropic muon pion.
//			std::cout<<"detector: muB mu pi"<<std::endl;
			break;



		}
		break;
	default:
		std::cout<<"You shouldn't be able to get here."<<std::endl;
		return -1;
}

// Some sensibility checks


if(channel_flag == CHAN_ELECPI && mS < MPION+ME )
{
	std::cout<<"ERROR: electron pion channel can't work with sterile masses below "<<MPION+ME<<" GeV."<<std::endl;
	return -1;
}
if(channel_flag == CHAN_MUONPI && mS < MMU+MPION)
{
	std::cout<<"ERROR: muon pion channel can't work with sterile masses below "<<MMU+MPION<<" GeV."<<std::endl;
	return -1;
}
if(channel_flag == CHAN_NUPI0 && mS < MPI0)
{
	std::cout<<"ERROR: nu pi0 channel can't work with sterile masses below "<<MPI0<<" GeV."<<std::endl;
	return -1;
}
if(channel_flag == CHAN_ELECPOSI && mS < 2*ME)
{
	std::cout<<"ERROR: ee channel can't work with sterile masses below "<<2*ME<<" GeV."<<std::endl;
	return -1;
}
if(channel_flag == CHAN_MUMU && mS < 2*MMU)
{
	std::cout<<"ERROR: mu mu channel can't work with sterile masses below "<<2*MMU<<" GeV."<<std::endl;
	return -1;
}
if(channel_flag == CHAN_MUE && mS < ME+MMU)
{
	std::cout<<"ERROR: mu e channel can't work with sterile masses below "<<MMU+ME<<" GeV."<<std::endl;
	return -1;
}

/***************************
 *	Main program Flow
 * ***********************/

if(matrix_flag == 0)
{
	switch(channel_flag) //this switch is currently pointless, but allows for channel dependent processing.
	{
		case CHAN_ELECPOSI:	
		case CHAN_ELECPI:
		case CHAN_MUONPI:	
			output_distributions(r, DETECTOR, CHAN, mS,hepevt, flux_name, NUMEVENTS);
			break;
		case CHAN_NUPI0:
		case CHAN_GAMMA:
		case CHAN_MUMU:
		case CHAN_MUE:	
			output_distributions(r, DETECTOR, CHAN, mS, hepevt, flux_name, NUMEVENTS);
			break;
		default:
			std::cout<<"Something bad has occured. Worry."<<std::endl;
			return -1;
	}
}


//Cleaning up.
gsl_rng_free(r);
delete CHAN;
delete DETECTOR;

return 0;
}

/*****************************************************************************
 *
 *		What follows is the main bit of the code of actual producing events
 *
 *
 * ***************************************************************************/


int output_distributions(gsl_rng * r, detector * DETECTOR, twoIP_channel * CHAN, double mS, bool hepevt, std::string flux_name, int NUMEVENTS)
{


fluxfile flux(flux_name, mS);


//static double events[NUMEVENTS][2]; //define the storage for all the events, [0] = E_s, [1] = cos\theta_s

//getEvents(mS,DETECTOR->flag,events); // we are always loading up the events from the standard BNB flux. 

//We enter the main loop over events. For each one, computing the relevant
//observables.
//


static OBSERVABLES Obs; //This struct is contained in "decay.h"; it specifically gives variables for a two body event (e+,e-)
Obs.chan_identifier= CHAN->chan_identifier; // needs this so it knows what resolutions to use in smearing

int m; 
for(m=0;m<NUMEVENTS;m++) 
{

	//The flux files I have don't provide phi or cos angles for the steriles, so
	//we generate them here. Also set up your positions and timings
	
	double	cosS = 0.99999;
	double  phiS = 2.0*M_PI*gsl_rng_uniform(r);
	double  enS  = flux.get_event(r); 	

	double posx=2.33*100*gsl_rng_uniform(r);
	double posy=2.56*100*gsl_rng_uniform(r);
	double posz=10.37*100*gsl_rng_uniform(r);
	double time =0;



		initial_sterile nus(mS, enS, cosS, phiS);

//		nus.labframeP.print("nus.labframeP");

		//We call the appropriate functions from the channels.
				switch(CHAN->chan_identifier){
					case CHAN_MUMU:
						CHAN->decayfunctionMassive(nus,MMU,MMU,0.0);
						break;
					case CHAN_MUE:
						CHAN->decayfunctionMassive(nus,MMU,ME,0.0);
						break;	
					case CHAN_ELECPI:
					case CHAN_MUONPI:
					case CHAN_NUPI0:
					case CHAN_ELECPOSI:
					case CHAN_GAMMA:
						CHAN->decayfunction(nus);
						break;
			}
		

		CHAN->observables(&Obs, r);

		// The following sterile observables can't be assigned at the channel
		// level anymore... in some sense they are inputs, not properties of the
		// outgoing event, so I'm not sure I think this is a problem.
		Obs.E_sterile = nus.energy; 	
		Obs.Th_sterile = nus.costhS;

	

		if(DETECTOR->accept(&Obs)==ACCEPTED)
		{


		if(hepevt == false){
			printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5f %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", Obs.E_sum, Obs.Th_sum, Obs.AngSep, Obs.E_sterile, Obs.Th_sterile, Obs.E_high, Obs.Th_high, Obs.E_low, Obs.Th_low, Obs.FS_AngSep,Obs.P_high[0],Obs.P_high[1],Obs.P_high[2],Obs.P_low[0],Obs.P_low[1],Obs.P_low[2],  Obs.E_sum_smear, Obs.Th_sum_smear, Obs.AngSep_smear, Obs.E_sterile_smear, Obs.Th_sterile_smear, Obs.E_high_smear, Obs.Th_high_smear, Obs.E_low_smear, Obs.Th_low_smear, Obs.FS_AngSep_smear, Obs.Minvar_smear);

			//for two body, E_high is ma, and E_low mb
		}
	       	if (hepevt)
		{
			
			switch(CHAN->chan_identifier){
				case CHAN_ELECPI:
				printf("%i 2\n",m);
				printf("1 11 0 0 0 0 %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",Obs.P_high[0],Obs.P_high[1],Obs.P_high[2],Obs.E_high,CHAN->model_params[0],posx,posy,posz,time);
				printf("1 211 0 0 0 0 %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",Obs.P_low[0],Obs.P_low[1],Obs.P_low[2],Obs.E_low,CHAN->model_params[1],posx,posy,posz,time);
				case CHAN_MUONPI:
				printf("%i 2\n",m);
				printf("1 13 0 0 0 0 %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",Obs.P_high[0],Obs.P_high[1],Obs.P_high[2],Obs.E_high,CHAN->model_params[0],posx,posy,posz,time);
				printf("1 211 0 0 0 0 %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",Obs.P_low[0],Obs.P_low[1],Obs.P_low[2],Obs.E_low,CHAN->model_params[1],posx,posy,posz,time);
				break;
				case CHAN_NUPI0:
				case CHAN_ELECPOSI:
				case CHAN_GAMMA:
				default:
					std::cout<<"ERROR: HEPevt not set up for this channel yet. Only e-pi and mu-pi."<<std::endl;
					return 0;
			}
		}
				// id# numpartic
				// track? PDGid mother1 mother2 daght1 daug2 px py pz En Mn Posx Posy Posy time


		}
	
}


return 0;

}

