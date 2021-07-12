#include "basic_definitions.h"

/****function to calculate timestep****/

void Tstep(int &imax, int &ib2, double *&x, double *&a, double *&vol, double **&cv, double *&p, double *&dt){
	
	double rho, u, cs, dx, sprad;

	for (int j=2; j<=ib2; j++){
		rho   = cv[0][j]/a[j];
		u     = cv[1][j]/cv[0][j];
		cs    = sqrt(Gamma*p[j]/rho);
		dx    = 0.5*(x[j+1]-x[j-1]);
		sprad = cs*sqrt(dx*dx+a[j]*a[j]) + fabs(u)*a[j];
		dt[j] = vol[j]/sprad;
	}

	dt[1]    = dt[2];
	dt[imax] = dt[ib2];

}
