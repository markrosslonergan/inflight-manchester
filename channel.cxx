#include "channel.h"

twoIP_channel::twoIP_channel(gsl_rng * g, std::vector<double> input_params)
{
	std::vector<double> p;
	p.push_back(0.0);	
	p.push_back(0.0);	
	p.push_back(0.0);	

	IP1.populate(1,p);
	IP2.populate(1,p);
	
	model_params = input_params;
	chan_identifier = model_params[2];	
	r = g;
}

int twoIP_channel::decayfunction(initial_sterile nuS)
{
	std::cout<<"You've somehow managed to call the decayfunction of the parent class (twoIP_channel). Don't do that."<<std::endl;
return 0;
}

int twoIP_channel::observables(OBSERVABLES * output, gsl_rng *g)
{
	//OBSERVABLES { double E_sum; double Th_sum; double AngSep; double E_sterile; double Th_sterile; double E_high; double Th_high; double E_low; double Th_low; double FS_AngSep; } OBSERVABLES;

	fourmomentum sum;
	std::vector<double> p1 = {IP1.p.at(0),IP1.p.at(1),IP1.p.at(2)};
	std::vector<double> p2 = {IP2.p.at(0),IP2.p.at(1),IP2.p.at(2)};

	bool IS_ORDER = false;

	output->P_low=p2;
	output->P_high=p1;

//		std::cout<<"test "<<(180.0/M_PI)*IP2.direction().at(0)<<" "<<acos(p2[2]/( pow(p2[0],2)+pow(p2[1],2)+pow(p2[2],2) ) )*180/3.14145<<
	
//			std::cout<<Obs.Th_low<<" "<<atan2(Obs.P_low[0],Obs.P_low[2])*180/3.14159<<" "<<atan2(Obs.P_low[1],Obs.P_low[2])*180/3.14159<<std::endl;
//			std::cout<<Obs.Th_low<<" "<<atan2(Obs.P_low[1],Obs.P_low[0])*180/3.14159<<" "<<atan2(Obs.P_low[2],Obs.P_low[0])*180/3.14159<<std::endl;
//			std::cout<<Obs.Th_low<<" "<<atan2(Obs.P_low[0],Obs.P_low[1])*180/3.14159<<" "<<atan2(Obs.P_low[1],Obs.P_low[2])*180/3.14159<<std::endl;
//			std::cout<<Obs.Th_low<<" "<<acos(Obs.P_low[2]/( pow(Obs.P_low[0],2)+pow(Obs.P_low[1],2)+pow(Obs.P_low[2],2) ) )*180/3.14145<<std::endl;



	std::vector<double> sum_p;
	sum_p.push_back(IP1.p.at(0) + IP2.p.at(0));
	sum_p.push_back(IP1.p.at(1) + IP2.p.at(1));
	sum_p.push_back(IP1.p.at(2) + IP2.p.at(2));
	
	sum.populate(IP1.E + IP2.E, sum_p);

	output->E_sum = sum.E;	
	output->Th_sum =(180.0/M_PI)*sum.direction().at(0);	

	output->AngSep =(180.0/M_PI)*acos((IP1.p.at(0)*IP2.p.at(0) + IP1.p.at(1)*IP2.p.at(1) + IP1.p.at(2)*IP2.p.at(2))/(IP1.modp*IP2.modp)); 

//	IP1.print("IP1");
//	IP2.print("IP2");

//	if(IP1.p.at(2) > 0 && IP2.p.at(2) > 0)
//	{
//		output->FS_AngSep = (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
//	}
//	else 
//	{
//		output->FS_AngSep = 180 - (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
//	}

	if(IP1.p.at(2)*IP2.p.at(2) >= 0 && IP2.p.at(0)*IP2.p.at(0) >= 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(fabs(atan(IP1.p.at(0)/IP1.p.at(2))) - fabs(atan(IP2.p.at(0)/IP2.p.at(2))));
	}
	else if(IP1.p.at(2)*IP2.p.at(2) >= 0 && IP2.p.at(0)*IP2.p.at(0) < 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) + atan(IP2.p.at(0)/IP2.p.at(2)));
	}
	else if(IP1.p.at(2)*IP2.p.at(2) < 0 && IP2.p.at(0)*IP2.p.at(0) >= 0)
	{
		output->FS_AngSep = 180.0- (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) + atan(IP2.p.at(0)/IP2.p.at(2)));
	}
	else if(IP1.p.at(2)*IP2.p.at(2) < 0 && IP2.p.at(0)*IP2.p.at(0) < 0)
	{
		output->FS_AngSep = 180.0- (180.0/M_PI)*fabs(atan(IP1.p.at(0)/IP1.p.at(2)) - atan(IP2.p.at(0)/IP2.p.at(2)));
	}


	output->E_high = IP1.E;	
	output->Th_high = (180.0/M_PI)*IP1.direction().at(0);	
	output->E_low = IP2.E;	
	output->Th_low = (180.0/M_PI)*IP2.direction().at(0);	

	double temp = 0.0;

	//Only want to do this in e+e- scenario, NOT mu pi scenario! hmm
	if(IS_ORDER && output->E_high < output->E_low)
	{ 	
		temp = output->E_low; 
		output->E_low = output->E_high; 
		output->E_high = temp; 
		
		temp = output->Th_low;
		output->Th_low = output->Th_high;
		output->Th_high = temp;
	}

	double mlow = sqrt(pow(IP2.E,2)-p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]);
	double mhigh = sqrt(pow(IP1.E,2)-p1[0]*p1[0]-p1[1]*p1[1]-p1[2]*p1[2]);

	output->Minvar = sqrt(mlow*mlow+mhigh*mhigh+2.0*(IP2.E*IP1.E-p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]));
	// Lets smear all of the above defined variables by preapproved gaussians



	double res_low = 0.0;
	double res_high = 0.0;
	double res_ang = 1.0;

	switch(output->chan_identifier)
	{
		case CHAN_ELECPOSI:
			res_low = 0.15;
			res_high = 0.15;
			break;
		case CHAN_ELECPI:
			res_low = 0.06;
			res_high = 0.15;
			break;
		case CHAN_MUONPI:
			res_low = 0.06;
			res_high = 0.06;
			break;
		case CHAN_NUPI0:
			res_low = 0.15;
			res_high = 0.15;
			break;
		case CHAN_GAMMA:
			res_low = 0.15;
			res_high = 0.15;
			break;
	}
	
	output->E_high_smear = smear_energy(output->E_high, mhigh, res_high, g);
	output->E_low_smear = smear_energy(output->E_low,  mlow, res_low, g);
	output->E_sum_smear = output->E_high_smear+output->E_low_smear;
	output->Th_high_smear = smear_angle(output->Th_high,res_ang,g);
	output->Th_low_smear = smear_angle(output->Th_low,res_ang,g);
	output->Th_sterile_smear = smear_angle(output->Th_sterile,res_ang,g);
	output->AngSep_smear = smear_angle(output->AngSep,res_ang,g);	
	output->FS_AngSep_smear = smear_angle(output->FS_AngSep,res_ang,g);	
	output->E_sterile_smear = output->E_sterile;
	output->Th_sum_smear = output->Th_sum_smear;

	output->Minvar_smear = sqrt(mlow*mlow+mhigh*mhigh+2.0*(output->E_low_smear*output->E_high-sqrt(pow(output->E_high_smear,2)-mhigh*mhigh)* sqrt(pow(output->E_low_smear,2)-mlow*mlow)*cos(output->AngSep_smear* 3.14159/180.0)  ) );


//	std::cout<<output->Minvar<<" "<<sqrt(mlow*mlow+mhigh*mhigh+2.0*(output->E_low_smear*output->E_high-sqrt(pow(output->E_high_smear,2)-mhigh*mhigh)* sqrt(pow(output->E_low_smear,2)-mlow*mlow)*cos(output->AngSep_smear* 3.14159/180.0)  ) )<<std::endl;

//	IP1.print("IP1");
//	IP2.print("IP2");

return 0;
}


/* ###############################
   
   Below here we have a section for each channel.

   ############################### */


//This is the nu_s \to \nu e+ e- channel (off-shell Zprime).

threebody::threebody(gsl_rng *g, std::vector<double> input) : twoIP_channel(g, input)
{
//	if(model_params.size() != 1)
//	{ 
//		std::cout<<"ERROR: threebody decay channel set up with too many parameters."<<std::endl; 
//	}
}

int threebody::decayfunction(initial_sterile nuS)
{
	double mZprime = model_params.at(0);
	double restFrameParams[] = {0.0,0.0,0.0};
	drawRestFrameDist(r,nuS.mass,mZprime,restFrameParams); //this will populate the doubles.
	computeLabFrameVariables(nuS.mass, nuS.energy, nuS.costhS, nuS.phiS, restFrameParams);
return 0;
}

int threebody::computeLabFrameVariables(double mS, double Es, double costhS, double phiS, double restFrameParams[3])
{
	double Enu = restFrameParams[0];
	double Th = acos(restFrameParams[1]);
	double Phi = restFrameParams[2];

//	std::cout<<"Enu: "<<Enu<<" cosThnu: "<<cos(Th)<<" PhiNu: "<<Phi<<std::endl;
//	std::cout<<"Es: "<<Es<<" cosThS: "<<costhS<<" phiS: "<<phiS<<std::endl;

	double me = 0.00;//511;// THIS CAUSES ERRORS! 5.11e-6; //GeV
	
	double Ee = (mS - Enu)/2.0;
	double Pe = sqrt(Ee*Ee-me*me);
	double beta = sqrt(1-mS*mS/(Es*Es));
	double gamma = 1.0/sqrt(1.0-beta*beta);

//	printf("Enu: %.5lf Ee: %.5lf Pe: %.5lf beta: %.5lf gamma: %.5lf mS: %.5lf\n", Enu, Ee, Pe, beta, gamma, mS); 

	double alpha = 2.0*acos( Enu/(2.0*Pe) );
	double theta_plus = M_PI - Th - alpha/2.0;
	double theta_minus = M_PI - Th + alpha/2.0;

	double Pplus_E = gamma*(Ee + beta*Pe*cos(theta_plus));
	double Pminus_E = gamma*(Ee + beta*Pe*cos(theta_minus));
	double Pplus_x = Pe*sin(theta_plus)*cos(Phi);
	double Pminus_x = Pe*sin(theta_minus)*cos(Phi);
	double Pplus_y = Pe*sin(theta_plus)*sin(Phi);
	double Pminus_y = Pe*sin(theta_minus)*sin(Phi);
	double Pplus_z = gamma*(Pe*cos(theta_plus) + beta*Ee);
	double Pminus_z = gamma*(Pe*cos(theta_minus) + beta*Ee);

	double Pee[] = {(Pplus_x + Pminus_x)/2.0, (Pplus_y + Pminus_y)/2.0, (Pplus_z + Pminus_z)/2.0};
	double Pplus[] = {Pplus_x, Pplus_y, Pplus_z};
	double Pminus[] = {Pminus_x, Pminus_y, Pminus_z};

//	std::cout<<"PplusNorm: "<<Pplus_E*Pplus_E - Pplus_x*Pplus_x - Pplus_y*Pplus_y - Pplus_z*Pplus_z<<std::endl;
//	std::cout<<"PmiusNorm: "<<Pminus_E*Pminus_E - Pminus_x*Pminus_x - Pminus_y*Pminus_y - Pminus_z*Pminus_z<<std::endl;

//	std::cout<<"Pp1: "<<Pplus[0]<<" Pp2: "<<Pplus[1]<<" Pp3: "<<Pplus[2]<<std::endl;
//	std::cout<<"Pm1: "<<Pminus[0]<<" Pm2: "<<Pminus[1]<<" Pm3: "<<Pminus[2]<<std::endl;

//	std::vector<double> Vec_pplus(Pplus, Pplus + sizeof(Pplus)/sizeof(double));
//	std::vector<double> Vec_pminus(Pminus, Pminus + sizeof(Pminus)/sizeof(double));

//	IP1.populate(Pplus_E, Vec_pplus);
//	IP1.print("pre-rot IP1");
//	IP2.populate(Pminus_E, Vec_pminus);
//	IP2.print("pre-rot IP2");

	rotor(acos(costhS),phiS,Pee);
	rotor(acos(costhS),phiS,Pplus);
	rotor(acos(costhS),phiS,Pminus);
//	std::cout<<"Rotated!..."<<std::endl;
//	std::cout<<"Pp1: "<<Pplus[0]<<" Pp2: "<<Pplus[1]<<" Pp3: "<<Pplus[2]<<std::endl;
//	std::cout<<"Pm1: "<<Pminus[0]<<" Pm2: "<<Pminus[1]<<" Pm3: "<<Pminus[2]<<std::endl;

	std::vector<double> Vec_pplus2(Pplus, Pplus + sizeof(Pplus)/sizeof(double));
	std::vector<double> Vec_pminus2(Pminus, Pminus + sizeof(Pminus)/sizeof(double));

	IP1.populate(Pplus_E, Vec_pplus2);
//	IP1.print("post-rot IP1");
	IP2.populate(Pminus_E, Vec_pminus2);
//	IP2.print("post-rot IP2");

//std::cout<<std::endl;

return 0;
}

double threebody::pdf_function(double x, double y, double mS, double mZprime, void * pointer)
{
	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1.0-mu_s*mu_s);

	double invnorm;

	if(alpha < 0.01)
	{
		double invnorm_perturb_0 = (1.0/(1.0+alpha))*(7.0/4.0 + 41.0/60.0*alpha); 

		double invnorm_perturb_rest = -0.18333*pow(alpha,2.0)+0.22857*pow(alpha,3.0)-0.23274*pow(alpha,4.0)+0.22421*pow(alpha,5.0)-0.21190*pow(alpha,6.0)+0.19899*pow(alpha,7.0)-0.18662*pow(alpha,8.0)+0.17517*pow(alpha,9.0)-0.16475*pow(alpha,10.0)+0.15531*pow(alpha,11.0);
		invnorm = invnorm_perturb_0+invnorm_perturb_rest;
	}
	else 
	{		
		invnorm = (3.0/(2.0*pow(alpha,3.0)))*(2+alpha-3*pow(alpha,2.0))/(1.0+alpha) + (4.0*alpha*alpha-3.0)*(log(1+alpha)/log(exp(1.0)))/pow(alpha,4.0);
	}

	double ret = (1.0/invnorm)*x*(4-x*x)/((1.0+alpha*x)*(1.0+alpha*x));

//	printf("inv. norm: %.5lf\talpha: %.5lf\tret: %.5lf\n",invnorm,alpha,ret);

	if(ret<0){ ret = 0.0; }

return ret; 
}


int threebody::rotor(double theta, double phi, double vector[3])
{
	double x=vector[0];
	double y=vector[1];
	double z=vector[2];

	double rdotn1 = x*cos(phi) + y*sin(phi);
	double rdotn2 = z;
	double rdotn3 = -x*sin(phi) + y*cos(phi);

	vector[0]=(cos(theta)*rdotn1 + sin(theta)*rdotn2)*cos(phi) - rdotn3*sin(phi);
	vector[1]=(cos(theta)*rdotn1 + sin(theta)*rdotn2)*sin(phi) + rdotn3*cos(phi);
	vector[2]=-sin(theta)*rdotn1 + cos(theta)*rdotn2;

return 0;
}

int threebody::drawRestFrameDist(gsl_rng * r, double mS, double mZprime, double output[3])
{

	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1-mu_s*mu_s);

//	double PDF_MAX = 3.0/2.0; //pdf_function_test
	double PDF_MAX = 1.8; //pdf_function

	double x = gsl_rng_uniform(r);
	double y = -1.0 + 2.0*gsl_rng_uniform(r);
	double phi = 2*M_PI*gsl_rng_uniform(r);
	double z = (PDF_MAX+0.01)*gsl_rng_uniform(r);

	while(threebody::pdf_function(x,y,mS,mZprime,NULL)<=z)
	{
//		printf("I tried!\n");
		x = gsl_rng_uniform(r);
		y = -1.0 + 2.0*gsl_rng_uniform(r);
		z = (PDF_MAX+0.01)*gsl_rng_uniform(r);
	}

//	printf("I succeeded!\n");

	//printf("%.5lf %.5lf %.5lf %.5lf\n",x,y,z,pdf(x,y,mS,mZprime,NULL));
	struct threebody::PDF_CHOICE choice;
	choice.Enu = mS*x/2.0;
	choice.cosThnu = y;
	choice.Phinu = phi;

	output[0]=choice.Enu;
	output[1]=choice.cosThnu;
	output[2]=choice.Phinu;

return 0;
}

//This is the nu_s \to \nu Zprime \to \nu e+ e- channel (on-shell Zprime).

Zprimeresonance::Zprimeresonance(gsl_rng *g, std::vector<double> input) : twoIP_channel(g, input)
{
//	if(model_params.size() != 1)
//	{
//	 	std::cout<<"ERROR: Zprime resonance channel set up with too many parameters."<<std::endl;
//	}
}

double Zprimeresonance::fourvec_costheta(double FOURVEC[4])
{
return FOURVEC[3]/sqrt(FOURVEC[1]*FOURVEC[1]+FOURVEC[2]*FOURVEC[2]+FOURVEC[3]*FOURVEC[3]);
}

double Zprimeresonance::fourvec_cosphi(double FOURVEC[4])
{
	double P = sqrt(FOURVEC[1]*FOURVEC[1]+FOURVEC[2]*FOURVEC[2]+FOURVEC[3]*FOURVEC[3]);
	double cosPhi = FOURVEC[1]/(sqrt(1-fourvec_costheta(FOURVEC)*fourvec_costheta(FOURVEC))*P);
return cosPhi;
}

double Zprimeresonance::rot_boost(double costheta, double phi, double gamma, double FOURVEC[4])
{
	double sintheta = sqrt(1-costheta*costheta);

	double beta = sqrt(1.0-1.0/(gamma*gamma));
	double temp[4];
	temp[0]=FOURVEC[0];
	temp[1]=FOURVEC[1];
	temp[2]=FOURVEC[2];
	temp[3]=FOURVEC[3];

	FOURVEC[0] = gamma*temp[0] + gamma*beta*temp[3];
	FOURVEC[1] = gamma*beta*cos(phi)*sintheta*temp[0] + cos(phi)*costheta*temp[1] - sin(phi)*temp[2] + gamma*cos(phi)*sintheta*temp[3];
	FOURVEC[2] = gamma*beta*sin(phi)*sintheta*temp[0] + sin(phi)*costheta*temp[1] + cos(phi)*temp[2] + gamma*sin(phi)*sintheta*temp[3];
	FOURVEC[3] = gamma*beta*costheta*temp[0] - sintheta*temp[1] + gamma*costheta*temp[3];

return 0.0;
}



int Zprimeresonance::decayfunction(initial_sterile nuS)
{

double mZprime = model_params.at(0);
double Es = nuS.energy;
double costhS = nuS.costhS;
double phiS = nuS.phiS;
double mS = nuS.mass;

//

double Z_E_srf = (mS*mS+mZprime*mZprime)/(2.0*mS);
double Z_P_srf = sqrt(Z_E_srf*Z_E_srf-mZprime*mZprime);

double Z_phi_srf = 0.0;
double Z_costheta_srf = 1.0;

double S_phi_lf = 0.0;
double S_costheta_lf = 1.0;
double S_E_lf = Es;
double S_gamma = Es/mS;

double Z_FOURVEC[] = {0.0,0.0,0.0,0.0};
double EPLUS_FOURVEC[] = {0.0,0.0,0.0,0.0};
double EMINUS_FOURVEC[] = {0.0,0.0,0.0,0.0};
double SUM_FOURVEC[] = {0.0,0.0,0.0,0.0};

double Z_gamma = 1.0;
double cosalpha, sinalpha;
double cosbeta, sinbeta, beta;
double temp = 0.0;

	// Angles of the Z in the sterile rest frame (srf) are evenly distributed on the sphere.
	Z_phi_srf = 2.0*M_PI*gsl_rng_uniform(r);
	Z_costheta_srf = 2.0*gsl_rng_uniform(r) -1.0;
	
	// The labframe phi and costheta for sterile (S_phi_lf, S_costheta_lf) are taken from input.
	S_phi_lf = phiS;
	S_costheta_lf = costhS;	
	S_E_lf = Es;
	S_gamma = S_E_lf/mS;

	// We define the Z fourvector in the sterile restframe.
	Z_FOURVEC[0] = Z_E_srf;
	Z_FOURVEC[1] = Z_P_srf*sqrt(1.0-Z_costheta_srf*Z_costheta_srf)*cos(Z_phi_srf);
	Z_FOURVEC[2] = Z_P_srf*sqrt(1.0-Z_costheta_srf*Z_costheta_srf)*sin(Z_phi_srf);
	Z_FOURVEC[3] = Z_P_srf*Z_costheta_srf;

	//We boost and rotate to move the Z fourvector into the lab frame.
	rot_boost(S_costheta_lf,S_phi_lf,S_gamma,Z_FOURVEC);
	
//	printf("%.5lf %.5lf %.5lf %.5lf %.5lf\n", Z_FOURVEC[0],  fourvec_costheta(Z_FOURVEC), fourvec_cosphi(Z_FOURVEC), Z_costheta_srf, cos(Z_phi_srf)); 

	Z_gamma = Z_FOURVEC[0]/mZprime;
	
	cosalpha = 2.0*gsl_rng_uniform(r) - 1.0;
	beta = 2.0*M_PI*gsl_rng_uniform(r);
	cosbeta = cos(beta);
	sinalpha = sqrt(1.0-cosalpha*cosalpha);
	sinbeta = sqrt(1.0-cosbeta*cosbeta);

	EPLUS_FOURVEC[0] = mZprime/2.0;
	EPLUS_FOURVEC[1] = (mZprime/2.0)*sinalpha*cosbeta;
	EPLUS_FOURVEC[2] = (mZprime/2.0)*sinalpha*sinbeta;
	EPLUS_FOURVEC[3] = (mZprime/2.0)*cosalpha;

	EMINUS_FOURVEC[0] = mZprime/2.0;
	EMINUS_FOURVEC[1] = -(mZprime/2.0)*sinalpha*cosbeta;
	EMINUS_FOURVEC[2] = -(mZprime/2.0)*sinalpha*sinbeta;
	EMINUS_FOURVEC[3] = -(mZprime/2.0)*cosalpha;

	rot_boost(fourvec_costheta(Z_FOURVEC), acos(fourvec_cosphi(Z_FOURVEC)),Z_gamma, EPLUS_FOURVEC); 
	rot_boost(fourvec_costheta(Z_FOURVEC), acos(fourvec_cosphi(Z_FOURVEC)),Z_gamma, EMINUS_FOURVEC); 

	SUM_FOURVEC[0] = EPLUS_FOURVEC[0] + EMINUS_FOURVEC[0], 
	SUM_FOURVEC[1] = EPLUS_FOURVEC[1] + EMINUS_FOURVEC[1];
	SUM_FOURVEC[2] = EPLUS_FOURVEC[2] + EMINUS_FOURVEC[2];
	SUM_FOURVEC[3] = EPLUS_FOURVEC[3] + EMINUS_FOURVEC[3];

	std::vector<double> eplus_p;
	eplus_p.push_back(EPLUS_FOURVEC[1]);
	eplus_p.push_back(EPLUS_FOURVEC[2]);
	eplus_p.push_back(EPLUS_FOURVEC[3]);

	std::vector<double> eminus_p;
	eminus_p.push_back(EMINUS_FOURVEC[1]);
	eminus_p.push_back(EMINUS_FOURVEC[2]);
	eminus_p.push_back(EMINUS_FOURVEC[3]);

	IP1.populate(EPLUS_FOURVEC[0], eplus_p);
	IP2.populate(EMINUS_FOURVEC[0], eminus_p);

return 0;
}

//This is a generic nu_s \to two body channel (isotropic in rest-frame decay)

twobody::twobody(gsl_rng * r, std::vector<double> input) : twoIP_channel(r, input)
{

//	if(input.size() != 2) 
//	{ 
//		std::cout<<"ERROR: twobody channel needs 2 inputs (both of the final particle masses)!"<<std::endl; 
//	}	

}

int twobody::decayfunction(initial_sterile nus)
{	
	//Pull the final state particle masses from the model_params input provided on initialization
	double Ma = model_params.at(0); 
	double Mb = model_params.at(1); 

	double Ea_rf = (nus.mass*nus.mass + Ma*Ma - Mb*Mb)/(2.0*nus.mass);
	double Eb_rf = (nus.mass*nus.mass + Mb*Mb - Ma*Ma)/(2.0*nus.mass);

//	std::cout<<"Ea_rf = "<<Ea_rf<<std::endl;
//	std::cout<<"Eb_rf = "<<Eb_rf<<std::endl;
	
	//both decay products have a common *magnitude* of momentum
	double P_rf = sqrt(Ea_rf*Ea_rf - Ma*Ma);

	//We assumed isotropy in the restframe... so we can generate angles flatly on the 2-sphere.
	double costheta_rf = 2.0*gsl_rng_uniform(r) - 1.0;
	double phi_rf  = 2.0*M_PI*gsl_rng_uniform(r);

	double sintheta_rf = sqrt(1-costheta_rf*costheta_rf);

	std::vector<double> momentum_a_rf;
	momentum_a_rf.push_back(P_rf*sintheta_rf*cos(phi_rf));
	momentum_a_rf.push_back(P_rf*sintheta_rf*sin(phi_rf));
	momentum_a_rf.push_back(P_rf*costheta_rf);

	//By momentum conservation this is just the negative of momentum_b;
	std::vector<double> momentum_b_rf;
	momentum_b_rf.push_back(-P_rf*sintheta_rf*cos(phi_rf));
	momentum_b_rf.push_back(-P_rf*sintheta_rf*sin(phi_rf));
	momentum_b_rf.push_back(-P_rf*costheta_rf);

	//Set up the decay product four momenta in the sterile rest frame.
	IP1.populate(Ea_rf, momentum_a_rf);
	IP2.populate(Eb_rf, momentum_b_rf);

//	IP1.print("pre-rot IP1");
//	IP2.print("pre-rot IP2");

	//We now need to boost and rotate so that the z axis lies along the sterile direction at an appropriate boost.
	IP1.rot_boost_from_parent(&nus.labframeP); 	
	IP2.rot_boost_from_parent(&nus.labframeP); 	

//	IP1.print("post-rot IP1");
//	IP2.print("post-rot IP2");

return 0;
}

