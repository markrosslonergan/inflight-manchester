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



#include "fourmomentum.h"


fourmomentum::fourmomentum(double energy, std::vector<double> momentum)
{
	populate(energy, momentum);
}

//This function allows the declaration of empty fourmometa to be filled by .populate(...) later. USE WITH CARE! Most of the functions would error with such a vector. 
fourmomentum::fourmomentum()
{
	E = 0.0;
	modp=0;
	p.push_back(0.0);
	p.push_back(0.0);
	p.push_back(0.0);
	mass = 0.0;
	
}

int fourmomentum::populate(double energy, std::vector<double> momentum)
{

	E = energy;
	p = momentum;

	//Just check that the input is OK.
	if(p.size() != 3) { std::cout<<"ERROR: 3-momentum wrong size."<<std::endl; }

	//Define the total 3-momentum and the full Minkowski norm (and check it's an on-shell 4-momentum).
	modp = sqrt(p.at(0)*p.at(0) + p.at(1)*p.at(1) + p.at(2)*p.at(2));

	mass = E*E - modp*modp; //Using mass as a temporary variable here. True value created a few lines down.
	if(fabs(mass) < 1e-12){ mass = 0.0; }

	if(mass < 0.0 ){ std::cout<<"ERROR: 4-vector is spacelike. This isn't what we had agreed on!"<<std::endl; }
	else{ mass = sqrt(mass); }


}

int fourmomentum::print(std::string name)
{

std::cout<<"Fourvector '"<<name<<"'"<<" = ("<<E<<", "<<p.at(0)<<", "<<p.at(1)<<", "<<p.at(2)<<"),\t"<<"[Inv. Mass^2: "<<mass<<", Norm of 3-momentum: "<<modp<<"]"<<std::endl; 

return 0;
}

std::vector<double> fourmomentum::direction()
{
	std::vector<double> temp;

	if(modp==0)
	{
		std::cout<<"Cannot compute direction: 3-momentum vanishes."<<std::endl;
		temp.push_back(0.0);
		temp.push_back(0.0);
	}
	else{
		temp.push_back(acos(p.at(2)/modp)); //acos returns 0 to pi.
		double phi = fabs(atan(p.at(1)/p.at(0)));

		if(p.at(0)>=0.0 && p.at(1) >= 0.0)
		{
			temp.push_back(phi);
		}
		else if(p.at(0)>=0.0 && p.at(1) < 0.0)
		{
			temp.push_back(-phi);
		}
		else if(p.at(0)<0.0 && p.at(1) >= 0.0)
		{
			temp.push_back(M_PI-phi);
		}
		else if(p.at(0)<0.0 && p.at(1) < 0.0)
		{
			temp.push_back(-M_PI+phi);
		}
		

	}

return temp;	
}

double fourmomentum::gamma()
{
double temp;
	if(mass==0)
	{ 
		std::cout<<"ERROR: Trying to compute gamma factor for NULL four momentum."<<std::endl; 
		temp = 1e-5;
	}	
	else
	{
		temp = E/mass;
	}

return temp;
}

int fourmomentum::rot_boost_from_parent(fourmomentum * parentP)
{
	double costheta = cos(parentP->direction().at(0));
	double sintheta = sin(parentP->direction().at(0));
	double phi = parentP->direction().at(1);

	double gamma = parentP->gamma();
	double beta = sqrt(1.0-1.0/(gamma*gamma));

	double temp[4];
	temp[0]=E;
	temp[1]=p.at(0);
	temp[2]=p.at(1);
	temp[3]=p.at(2);

	double new_E;
	std::vector<double> new_p;

	new_E = gamma*temp[0] + gamma*beta*temp[3];
	new_p.push_back( gamma*beta*cos(phi)*sintheta*temp[0] + cos(phi)*costheta*temp[1] - sin(phi)*temp[2] + gamma*cos(phi)*sintheta*temp[3] );
	new_p.push_back( gamma*beta*sin(phi)*sintheta*temp[0] + sin(phi)*costheta*temp[1] + cos(phi)*temp[2] + gamma*sin(phi)*sintheta*temp[3] );
	new_p.push_back( gamma*beta*costheta*temp[0] - sintheta*temp[1] + gamma*costheta*temp[3] );

	fourmomentum::populate(new_E, new_p);

return 0;
}



