#include "basic_definitions.h"

/****forward declaration of entropy correction function****/
double Entropy_corr(double &z, double &d);

void Flux_roe(int &imax, int &ib2, double *&a, double **&ls, double **&rs){
	double **f;
	double *fcav, *fdiss;
	double ggm1, si, rl, ul, pl, hl, qrl, rr, ur, pr, hr, qrr, rav, dd, dd1, uav, hav, q2a, c2a, cav, du, delta,
	       h1, h2, h3, h5, eabs1, eabs2, eabs3;

	fcav = new double[5]; fdiss = new double[5];
	Allocate_2D(f, noconv, idim);

	ggm1=Gamma/gam1;
	
	
	for(int j=1; j<=ib2; j++){
		si = 0.5*(a[j+1]+a[j]);

/**** average of left and right state****/
		rl  = ls[0][j];
		ul  = ls[1][j];
		pl  = ls[2][j];
		hl  = ggm1*pl/rl + 0.5*ul*ul;
		qrl = ul*rl;

		rr  = rs[0][j]; 
		ur  = rs[1][j];
		pr  = rs[2][j];
		hr  = ggm1*pr/rr + 0.5*ur*ur;
		qrr = ur*rr;
		
		fcav[0] = qrl + qrr;
		fcav[1] = qrl*ul + qrr*ur + pl + pr;
		fcav[2] = qrl*hl + qrr+hr;

/**** dissipative flux****/

		rav    = sqrt(rl*rr);
		dd     = rav/rl;
		dd1    = 1.0/(1.0+dd);
		uav    = (ul+dd*ur)*dd1;
		hav    = (hl+dd*hr)*dd1;
		q2a    = 0.5*uav*uav;
		c2a    = gam1*(hav-q2a);
		cav    = sqrt(c2a);
		du     = ur-ul;

		h1     = fabs(uav-cav);
		h2     = fabs(uav);
		h3     = fabs(uav+cav);
		delta  = epsentr*cav;
		
		eabs1  = Entropy_corr(h1, delta);
		eabs2  = Entropy_corr(h2, delta);
		eabs3  = Entropy_corr(h3, delta);
		
		h1     = rav*cav*du;
		h2     = eabs1*(pr-pl-h1)/(2.0*c2a);
		h3     = eabs2*(rr-rl-(pr-pl)/c2a);
		h5     = eabs3*(pr-pl+h1)/(2.0*c2a);
		
		fdiss[0] = h2+h3+h5;
		fdiss[1] = h2*(uav-cav)+h3*uav+h5*(uav+cav);
		fdiss[2] = h2*(hav-cav*uav)+h3*q2a+h5*(hav+cav*uav);

/****total flux at i+1/2****/
		f[0][j] = 0.5*(fcav[0]-fdiss[0])*si;
		f[1][j] = 0.5*(fcav[1]-fdiss[1])*si;
		f[2][j] = 0.5*(fcav[2]-fdiss[2])*si;
	}
/****summ of fluxes = RHS ****/
	
	for(int j=2; j<=ib2; j++){
		rhs[0][j] = f[0][j]-f[0][j-1];	
		rhs[1][j] = f[1][j]-f[1][j-1];
		rhs[2][j] = f[2][j]-f[2][j-1];
	}
}


/****entropy correction function to avoid sonic point and other problems****/
double Entropy_corr(double &z, double &d){
	if(z>d) return z;
	else return 0.5*(z*z+d*d)/d;
}


