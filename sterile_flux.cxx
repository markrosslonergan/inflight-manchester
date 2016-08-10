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

#include "sterile_flux.h"

initial_sterile::initial_sterile(double M, double E, double in_costhS, double in_phiS)
{
	mass = M;
	energy = E;
	costhS = in_costhS; 
	phiS = in_phiS;
	double totP =sqrt(energy*energy-M*M);	

	double temp[] = {totP*sqrt(1.0-in_costhS*in_costhS)*cos(in_phiS), totP*sqrt(1.0-in_costhS*in_costhS)*sin(in_phiS), totP*in_costhS };

	std::vector<double> momentum(temp, temp + sizeof(temp)/sizeof(double));

	labframeP.populate(energy, momentum);
}

double getEvents(double mS, int detector_flag, double events[][2])
{
	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;

	int n = 1;
	int m = 0;
	char s[100];
	char filename[500] = "MC_flux/\0";

	if(detector_flag == DET_NOCUTS)
	{
		sprintf(s,"%.4lf_MUBOONE.dat", mS);
	}
	else if(detector_flag == DET_SBND)
	{
		sprintf(s,"%.4lf_SBND.dat", mS);
	}
	else if(detector_flag == DET_MUBOONE)
	{
		sprintf(s,"%.4lf_MUBOONE.dat", mS);
	}
	else if(detector_flag == DET_ICARUS)
	{
		sprintf(s,"%.4lf_ICARUS.dat", mS);
	}
	else 
	{	std::cout<<"Look I'm not smart enough for this. (detector_flag = "<<detector_flag<<")"<<std::endl;
		return 0;
	}

	strcat(filename,s);
//	printf("Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS. sterileflux.cxx \n");
		exit(1);
	}
    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		pch = strtok(buf,"\t");
		n=0;
 		while (pch != NULL)
		{
			if(n==0){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][0] = strtof(pch,NULL);	 
				}
			if(n==1){	//printf("%.7g\n", strtof(pch,NULL));
					events[m][1] = strtof(pch,NULL);	 
				}
//		std::cout<<"n = "<<n<<std::endl;
			pch = strtok(NULL,"\t");
		n++;	
		}
		m++;
	}
	fclose(ptr_file);

//printf("Total lines: %d\n", m);

return 0;
}


fluxfile::fluxfile(std::string name, double mSin)
{

mass = mSin;

filename = name;

fmax = 0.0;

if(filename != "none"){

	int k = 0;
	std::string strE;
	std::string strW;
	std::ifstream myfile (filename);
	if (!myfile.is_open())
	{
		std::cout<<"#ERROR: flux::flux @flux.c, passed flux file does not exist"<<std::endl;
		exit(EXIT_FAILURE);
	}
		while(!myfile.eof()){
			
			myfile >> strE;
			myfile >> strW;
			double ien = atof(strE.c_str());
			double iflux = atof(strW.c_str());
			elist.push_back(ien);
			fluxlist.push_back(iflux);

			if(ien < mass && iflux != 0.0  ){
				std::cout<<"#ERROR: flux file containts events with energy less than the sterile rest mass! mass "<<mass<<" E "<<ien<<" flux "<<iflux<<std::endl;
				exit(EXIT_FAILURE);
			}

			if(iflux > fmax){
				fmax=iflux;
			}



	   	//	std::cout<<"E: "<<elist.back()<<" W: "<<fluxlist.back()<<std::endl;
		}	

		elist.pop_back();
		fluxlist.pop_back();



	myfile.close();

	if(elist.size()!=fluxlist.size()){
		std::cout<<"WHAT? elist != wlist in size. error @interpolate() in sterileflux.cxx. Your flux file has uneven columns"<<std::endl;
		exit(EXIT_FAILURE);
	}

//	for(int i =0; i<elist.size(); i++){
//		std::cout<<"elist: "<<i<<" "<<elist[i]<<" "<<fluxlist[i]<<std::endl;
//	}		

	}
}


double fluxfile::getFlux(double Evin)
{

	double ans = interpolate(elist, fluxlist, Evin);

	if(Evin < mass){ //we already know there is no flux below mass, checked when reading in. But linear interpolation doesnt know this. enforce here
		ans = 0.0;
	}


return ans;
}



double interpolate(std::vector<double > elist, std::vector<double> wlist, double Ein)
{
	if(elist.size()!=wlist.size())
	{

			std::cout<<"#ERROR:internal interpolate failure, vector of different size"<<std::endl;
			exit(EXIT_FAILURE);

	}



	double e0 = 0.0;
	double e1 = 0.0;
	double w0 = 0.0;
	double w1 = 0.0;
	
	if(Ein > elist.back() || Ein < elist.front())
	{
		return 0.0;
	} 
	else if(elist.front()==Ein)
	{
		return wlist.front();
	}
	else if(elist.back()==Ein)
	{
		return wlist.back();
	}

	int i = 1;
	for(i=1; i<=elist.size();i++)
	{
		if(elist[i]==Ein)
		{
			return wlist[i];
		}

		if(elist[i] > Ein)
		{
			e0 = elist[i-1];
			w0 = wlist[i-1];

			e1 = elist[i];
			w1 = wlist[i];
			break;
		}
	}

	double ans = w0 +(w1-w0)*(Ein - e0)/(e1-e0);
		


	return ans; 
}


double fluxfile::get_event(gsl_rng *r)
{

	double height_draw = 0.0;
	double norm_flux = 0.0;
	double E = 0.0;
	double current = 0;
	double E_max=5.0;
	double E_min=1e-5;

	int num = 0;

	E = (E_max-E_min)*gsl_rng_uniform(r)+E_min;

	height_draw = gsl_rng_uniform(r);
	norm_flux = getFlux(E)/fmax;

	while(height_draw > norm_flux)
	{
		num++;
		E = (E_max-E_min)*gsl_rng_uniform(r)+E_min;
		height_draw = gsl_rng_uniform(r);
		norm_flux = getFlux(E)/fmax;
	//	if(num > 1){std::cout<<"calling a lot "<<num<<"  E: "<<E<<" fmax: "<<fmax<<"  getFlux(E): "<<getFlux(E)<<" "<<norm_flux<<std::endl;}
	}



return E;
}







