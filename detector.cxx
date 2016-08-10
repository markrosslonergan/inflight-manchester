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

#include "detector.h"


double smear_energy(double En,double mn, double Percen, gsl_rng * r){
		double ans = 0;
		while(ans <=mn){
		ans =gsl_ran_gaussian ( r,Percen*En/sqrt(En))+En;
		}
		return ans;
	} 

double smear_angle(double th, double Percen,gsl_rng * r){
		double ans = 0;
		ans = gsl_ran_gaussian(r,Percen)+th;
		
		return fabs(ans);
	} 



//This is just to make sure that we have set up the detector object correctly.
int detector::accept(OBSERVABLES * Obs)
{
	std::cout<<"I'm the accept function of a ghost detector!"<<std::endl;

return -1;
}

//the pseudo-detector accepts all. It casts no judgement.
int nocuts::accept(OBSERVABLES * Obs)
{
	return ACCEPTED;
}


muBooNE::muBooNE()
{
	flag = DET_MUBOONE;
	std::cout<<"inside class det flag = "<<flag<<std::endl;

	Energy_threshold = 0.1; 	// In GeV
	AngSep_threshold = 30.0; 	// In Degrees
	Energy_ratio_threshold = 0.1; 	// Percentage.
}

int muBooNE::accept(OBSERVABLES * Obs)
{
    if(	Obs->E_sum < Energy_threshold || 
	Obs->FS_AngSep < AngSep_threshold ||
	Obs->E_low/Obs->E_high < Energy_ratio_threshold ) 
    { 
	return REJECTED; 
    }
    else
    {
	return ACCEPTED;
    }

}

//These below should be filled out once we know the acceptance criteria for
//muBooNE in all channels.
muBooNE_ee::muBooNE_ee()
{
	flag = DET_MUBOONE;

	AngSep_threshold_low = 2; //Below this electorn positron pairs are seen as a single track.
	AngSep_threshold_high = 40; //background cut on events that are seperated by more than this
	Energy_sum_threshold = 0.1 ; //onyl accept energetic enough events

}

int muBooNE_ee::accept(OBSERVABLES * Obs)
{
	double Eenergy = 0.0;
	double THenergy = 0.0;
	if(Obs->E_low > Obs->E_high){
		Eenergy = Obs->E_low_smear;
		THenergy = Obs->Th_low_smear;
	}else{
		Eenergy = Obs->E_high_smear;
		THenergy = Obs->Th_high_smear;
	}	


    if( THenergy > 20  || Obs->AngSep_smear < 2 || Eenergy < 0.1 ) 
    { 
	return REJECTED; 
    }
    else
    {
	return ACCEPTED;
    }
}
//*********************** Gamma ***************************************
muBooNE_gamma::muBooNE_gamma()
{
	flag = DET_MUBOONE;

	 Energy_threshold = 0.3; //Below this events aren't registered.
	 AngSep_threshold = 2.0; //Below this electorn positron pairs are seen as a single track.
	 Ang_threshold = 5.0; //

}

int muBooNE_gamma::accept(OBSERVABLES * Obs)
{	
    if(	Obs->AngSep > 2  ) 
    { 
	return REJECTED; 
    }
    else
    {
	return ACCEPTED;
    }
}

//######
//*********************** mu pi ***************************************
muBooNE_mupi::muBooNE_mupi()
{
	flag = DET_MUBOONE;

	 Energy_sum_threshold = 1.0; //Below this events aren't registered.
	 AngSep_threshold = 5.0; //Below this electorn positron pairs are seen as a single track.

}

int muBooNE_mupi::accept(OBSERVABLES * Obs)
{

    if(	Obs->AngSep_smear < 40.0 || Obs->Th_high_smear > 80 || Obs->Th_low_smear > 80 ) // || Obs->Th_high_smear > Obs->AngSep_smear ||  Obs->Th_low_smear > Obs->AngSep_smear  ) 
    { 
	return REJECTED; 
    }
    else
    {
	return ACCEPTED;
    }
}

//######
//######
//*********************** epi ***************************************
muBooNE_epi::muBooNE_epi()
{
	flag = DET_MUBOONE;

	 Energy_sum_threshold = 1.0; //Below this events aren't registered.
	 AngSep_threshold = 5.0; //Below this electorn positron pairs are seen as a single track.

}

int muBooNE_epi::accept(OBSERVABLES * Obs)
{	
    if(		Obs->AngSep_smear < 40 )//40.0 ) 
    { 
	return REJECTED; 
    }
    else
    {
	return ACCEPTED;
    }
}

//*********************** nupi0 ***************************************
muBooNE_nupi0::muBooNE_nupi0()
{
	flag = DET_MUBOONE;

	Energy_reco_threshold = 0.6; //Below this events aren't registered.
	Ang_reco_threshold =10; //
}

int muBooNE_nupi0::accept(OBSERVABLES * Obs)
{	
    if(	Obs->Th_low_smear > 25 /* not the ang_sep here, but angsep between photons || 
	Obs->AngSep_smear > 40*/ ) 
    { 
	return REJECTED; 
    }
    else
    {
	return ACCEPTED;
    }
}

//######
