#ifndef FOURVECTOR_H_
#define FOURVECTOR_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>


//This tacitly assumes that we are talking about on-shell 4-momenta... so
//timelike or null. I refuse to think about exceptions to this or how general
//the class could be/is.
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
