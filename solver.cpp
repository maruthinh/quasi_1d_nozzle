#include "basic_definitions.h"

/****where all the function required to solve equations are called****/

void Solver(ofstream &ofile_conv, int &iter, int &imax, int &ib2, int mxdum, double *&x, double *&a, double *&vol, double **&cv, double *&p, double **&cvold, double **&diss, double **&rhs, double *&dt, double *&dum){

	int irk, idim;	
	double fac, adtv, rrho, rhou, rhoe;

/****check the dimensions of dummy array (Dissp, Flux, LR_state, Irsmoo)****/

	idim = 4*imax;
	
	if(idim>mxdum){
		cout<<"maximum dimensions of work space exceeded in (solver.cpp)"<<endl;
	}
	
	for (int j=1; j<=imax; j++){
		cvold[0][j] = cv[0][j];
		cvold[1][j] = cv[1][j];
		cvold[2][j] = cv[2][j];
		diss[0][j]  = 0.0;
		diss[1][j]  = 0.0;
		diss[2][j]  = 0.0;
	}

/****calculate timestep****/

	Tstep(imax, ib2, x, a, vol, cv, p, dt);
	
/****loop over R-K stages****/

	for(irk=0; irk<nrk; irk++){

/**** flux difference splitting****/
		
		if(c[1]=='R'){
			LR_State_roe(imax, ib2, a, vol, cv, p);
			Flux_roe(imax, ib2, a, ls, rs);
		}
				
		else if(c[1]=='M'){
			LR_State_movers(imax, ib2, a, vol, cv, p);		
			Minmod_lim(imax, ib2, a, vol, cv);
			Flux_movers(imax, ib2, a, ls, rs, cvrs, cvls);		
		}
		
		else{
			cout<<"wrong input to the type of scheme to be used"<<endl;
			exit(0);
		}

/**** source term****/
		Srcterm(imax, ib2, a, p);

/**** residual timestep****/	
		fac = ark[irk]*cfl;
	
		
		for(int j=2;j<=ib2;j++){
			adtv      = fac*dt[j]/vol[j];
			rhs[0][j] = adtv*rhs[0][j];
			rhs[1][j] = adtv*rhs[1][j];
			rhs[2][j] = adtv*rhs[2][j];
		}

/****implicit residual smoothing****/

/****Update conserved variables and pressure****/
		
		for(int j=2;j<=ib2;j++){
			cv[0][j] = cvold[0][j] - rhs[0][j];
			cv[1][j] = cvold[1][j] - rhs[1][j];
			cv[2][j] = cvold[2][j] - rhs[2][j];

			rrho = a[j]/cv[0][j];
		        rhou = cv[1][j]/a[j];
		        rhoe = cv[2][j]/a[j];
		        p[j] = gam1*(rhoe-0.5*rhou*rhou*rrho);
			if (isnan(p[j])) {
				cout<<"\t"<<p[j]<<"\t"<<rrho<<"\t"<<rhou<<endl;
				exit(0);
			}
		}
		
		Bcond(imax, ib2, a, cv, p);

	}
	
	Conver(ofile_conv, iter, imax, ib2, a, cv, cvold, p);
}






















