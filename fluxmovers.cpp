#include "basic_definitions.h"

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min);

template <typename T>
T Max3(T a1, T a2, T a3);

template <typename T>
T Min3(T a1, T a2, T a3);



void Flux_movers(int &imax, int &ib2, double *&a, double **&ls, double **&rs, double **cvrs, double **cvls){
	double **f;
	double *fcav, *fdiss, *lambdar, *lambdal, *movers;
	double ggm1, si, rl, ul, pl, hl, qrl, rr, ur, pr, hr, qrr, lmax, lmin, er, el; 

	fcav = new double[3]; fdiss = new double[3]; lambdar = new double[3]; lambdal = new double[3];
	movers = new double[3];
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
		el  = pl/(rl*(Gamma-1.0)) + 0.5*ul*ul;

		rr  = rs[0][j]; 
		ur  = rs[1][j];
		pr  = rs[2][j];
		hr  = ggm1*pr/rr + 0.5*ur*ur;
		qrr = ur*rr;
		er  = pr/(rr*(Gamma-1.0)) + 0.5*ur*ur;
		
		fcav[0] = qrl + qrr;
		fcav[1] = qrl*ul + qrr*ur + pl + pr;
		fcav[2] = qrl*hl + qrr+hr;

/**** dissipative flux****/
		
		lambdar[0] = fabs(ur - sqrt(Gamma*pr/rr));
		lambdar[1] = fabs(ur);
		lambdar[2] = fabs(ur + sqrt(Gamma*pr/rr));

		lambdal[0] = fabs(ul - sqrt(Gamma*pl/rl));
		lambdal[1] = fabs(ul);
		lambdal[2] = fabs(ul + sqrt(Gamma*pl/rl));

		lmax = max(Max3(lambdar[0], lambdar[1], lambdar[1]), Max3(lambdal[0], lambdal[1], lambdal[2])); 
		lmin = max(Min3(lambdar[0], lambdar[1], lambdar[1]), Min3(lambdal[0], lambdal[1], lambdal[2])); 
	
		movers[0]=Movers(qrr, qrl, rr, rl, lmax, lmin);
		movers[1]=Movers(qrr*ur+pr, qrl*ul+pl, qrr, qrl, lmax, lmin);
		movers[2]=Movers(qrr*hr, qrl*hl, rr*er, rl*el, lmax, lmin);
		
		if(iorder==1){	
			fdiss[0] = (movers[0]+minmodlim[0][j]*(max(lambdar[1], lambdal[1])-movers[0]))*(rr-rl);
			fdiss[1] = (movers[1]+minmodlim[1][j]*(max(lambdar[1], lambdal[1])-movers[1]))*(qrr-qrl);
			fdiss[2] = (movers[2]+minmodlim[2][j]*(max(lambdar[1], lambdal[1])-movers[2]))*(rr*er-rl*el);
		}
		else if(iorder==2){	
			fdiss[0] = (movers[0]+minmodlim[0][j]*(max(lambdar[1], lambdal[1])-movers[0]))*(cvrs[0][j]-cvls[0][j]);
			fdiss[1] = (movers[1]+minmodlim[1][j]*(max(lambdar[1], lambdal[1])-movers[1]))*(cvrs[1][j]-cvls[1][j]);
			fdiss[2] = (movers[2]+minmodlim[2][j]*(max(lambdar[1], lambdal[1])-movers[2]))*(cvrs[2][j]-cvls[2][j]);
		}	
		else{
			cout<<"wrong input for order of accuracy"<<endl;
		}			
	

		/*fdiss[0] = movers[2]*(rr-rl);
		fdiss[1] = movers[2]*(qrr-qrl);
		fdiss[2] = movers[2]*(rr*er-rl*el);*/

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


/****function to compute movers dissipation****/
/*template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S;
	const double epsilon=1e-4;

	if (fabs(Fr-Fl)<=epsilon){
		if (fabs(Ur-Ul)<=epsilon) return L_min;
		else if (fabs(Ur-Ul)>epsilon) return fabs(((Fr-Fl)/(Ur-Ul)));		
	}
	
	else if(fabs(Fr-Fl)>epsilon){
		if (fabs(Ur-Ul)<=epsilon) S=L_min;
		else if (fabs(Ur-Ul)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
	}	
	
	if ((S)>=L_max) return (L_max);
	else if((S)<=L_min) return (L_min);
	else return (S);
}*/
/*
template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S, S1;
	const double epsilon=1e-4, delta=0.8;

	if (fabs(Fr-Fl)<=epsilon){
		if (fabs(Ur-Ul)<=epsilon) return L_min;
		else if (fabs(Ur-Ul)>epsilon) return fabs(((Fr-Fl)/(Ur-Ul)));		
	}
	
	else if(fabs(Fr-Fl)>epsilon){
		if (fabs(Ur-Ul)<=epsilon) S=L_min;
		else if (fabs(Ur-Ul)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
	}	
	
	if ((S)>=L_max) S1 = (L_max);
	else if((S)<=L_min) S1 = (L_min);
	else S1 = (S);
	
	if (fabs(S1)<delta) return (S1*S1+delta*delta)/(2.0*delta);
	else return S1;
}*/


/*template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S, delta=1e-5;

	if (Fr!=Fl and Ur!=Ul) S = fabs(((Fr-Fl)/(Ur-Ul)));
	else if (fabs(Fr-Fl)<delta) return 0.0;
	else return L_min;

	if ((S)>=L_max) return (L_max);
	else if((S)<=L_min) return (L_min);
	else return (S);
}*/

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S, delta=1e-2;

	if (fabs(Fr-Fl)<delta) return 0.0;
	else if (fabs(Ur-Ul)<delta) return L_min;
	else  S = fabs(((Fr-Fl)/(Ur-Ul)));

	if ((S)>=L_max) return (L_max);
	else if((S)<=L_min) return (L_min);
	else return (S);
}


/****to compute maximum value of 3 parameters****/
template <typename T>
T Max3(T a1, T a2, T a3){

	if (a1>=a2 && a1>=a3) return a1;
	else if (a2>=a3 && a2>=a1) return a2;
	else return a3;

}

/****to compute minimum value of 3 parameters****/
template <typename T>
T Min3(T a1, T a2, T a3){

	if (a1<=a2 && a1<=a3) return a1;
	else if (a2<=a3 && a2<=a1) return a2;
	else return a3;

}
