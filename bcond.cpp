#include "basic_definitions.h"

/*******************sets boundary conditions at the inlet and outlet of the nozzle************************/

void Bcond(int &imax, int &ib2, double *& a, double **& cv, double *& p){
	
	double u, cs2, c02, rinv, dis, cb, cc02, tb, pb, rhob, rho, cs, ub;

/***********************************************************************************************/
/* 			inlet subsonic							       */		
/*                                                                       		       */
/* 			speed of sound from Riemann invariant                                  */  
/* 			temperature from isentropic relation                                   */ 
/*		 	pressure from isentropic relation                                      */
/*		 	density from gas equation                                              */
/*		 	velocity from energy equation                                          */
/*                                                                                             */
/***********************************************************************************************/

	u    = cv[1][2]/cv[0][2];
	cs2  = Gamma*p[2]*a[2]/cv[0][2];
	c02  = cs2 + 0.5*gam1*u*u;
	rinv = u-2.0*sqrt(cs2)/gam1;
	dis  = gap1*c02/(gam1*rinv*rinv)-0.5*gam1;
	
	if(dis<0) {
		cout<<"dis in inflow boundary condition is less than zero so over writing"<<endl;
		dis=1e-20;				
	}
	
	cb   = -rinv*(gam1/gap1)*(1.0+sqrt(dis));
	cc02 = min(cb*cb/c02, 1.0);
	tb   = cc02*t01;
	pb   = p01*pow(tb/t01, Gamma/gam1);
	rhob = pb/(rgas*tb);
	ub   = sqrt(2.0*cpgas*(t01-tb));

	cv[0][1] = rhob*a[2];
	cv[1][1] = rhob*a[2]*ub;
	cv[2][1] = (pb/gam1+0.5*rhob*ub*ub)*a[2];
	p[1]     = pb;	

	//cout<<"to check boundary values"<<endl;
	//cout<<cv[1][1]<<"\t"<<cv[2][1]<<"\t"<<cv[3][1]<<"\t"<<p[1]<<endl;
// --- outlet - subsonic: 
//
//     - pressure = given back pressure
//     - density from characteristic b.c.
//     - velocity from characteristic b.c.
//
// --- outlet - supersonic:
//
//     - pressure extrapolated from the interior
//     - density extrapolated from the interior
//     - velocity extrapolated from the interior 

	rho = cv[0][ib2]/a[ib2];
	u   = cv[1][ib2]/cv[0][ib2];
	cs  = sqrt(Gamma*p[ib2]/rho);

	if (u>cs){
/****supersonic****/
		pb   = p[ib2];
		rhob = rho;
		ub   = u;
	}
	
	else{
/****subsonic****/
		pb   = p2;
		rhob = rho + (p2-p[ib2])/(cs*cs);
		ub   = u   - (p2-p[ib2])/(cs*rho);		
	}

	cv[0][imax] = rhob*a[ib2];
	cv[1][imax] = rhob*ub*a[ib2];
	cv[2][imax] = (pb/gam1 + 0.5*rhob*ub*ub)*a[ib2];
	p[imax]     = pb;
	
}


