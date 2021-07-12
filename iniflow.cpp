#include "basic_definitions.h"

/*Initializes the flow field using total pressure and temperature (inlet),
  and static pressure (outlet). Initializes reference values of the limiter
  (CUSP scheme or Roe's flux-difference splitting scheme)*/

/***************for expressions on total quantities see Modern Compressible flow by Anderson., page 59***************/


void Iniflow(int &imax, double *&a, double **&cv, double *&p){

	double temp, rho, mach, cs, u, mass, e;
	
	temp = t01*pow(p2/p01, gam1/Gamma);
	rho  = p2/(rgas*temp);
	mach = sqrt(2.0*((t01/temp)-1.0)/gam1);
	cs   = sqrt(Gamma*p2/rho);
	u    = cs*mach;
	mass = rho*u*a[2];
	e    = (cpgas-rgas)*t01;


/*********initialize flow field*************/

	for(int j=1; j<=imax; j++){
		cv[0][j] = rho*a[j];
		cv[1][j] = mass;
		cv[2][j] = rho*e*a[j];
 		p[j]     = p2;
	}
	
/**********limiter reference values************/
	
	volref = 1.0;
	rhoref = rho;
	uref   = u;
	pref   = p2;
}
		
