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


#ifndef FOURVECTOR_H_
#define FOURVECTOR_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>


//This tacitly assumes that we are talking about on-shell 4-momenta... so
//timelike or null
class fourmomentum {

public: 
	double mass;
	double E;
	double modp;
	std::vector<double> p;

	//For now I assume we define 4vectors at construction by their energies
	//and 3-momentum.  We could overload this for other combinations if
	//they look useful later on.
	fourmomentum(double energy, std::vector<double> momentum);
	fourmomentum();
	int populate(double energy, std::vector<double> momentum);
	int print(std::string name);
	std::vector<double> direction();
	double gamma();
	int rot_boost_from_parent(fourmomentum * parentP);
};


#endif
