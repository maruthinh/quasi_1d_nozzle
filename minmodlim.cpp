#include "basic_definitions.h"

template <typename T>
T MinMod(T Ujp2, T Ujp1, T Ujp);

template <typename T>
T Sign(T a);

/****computes left and right state using MUSCL appraoch****/

void Minmod_lim(int &imax, int &ib2, double *&a, double *&vol, double **&cv){
	
	double **du;

/****dynamic memory allocation for 2D & 1D pointers****/
	Allocate_2D(du, noconv, idim);

/****first difference of rho, u, p****/
	for(int j=1;j<=ib2;j++){
		du[0][j]    = cv[0][j+1]/a[j+1] - cv[0][j]/a[j];
		du[1][j]    = cv[1][j+1]/a[j+1] - cv[1][j]/a[j];
		du[2][j]    = cv[2][j+1]/a[j+1] - cv[2][j]/a[j];
	}
		du[0][0]    = du[0][1];
		du[1][0]    = du[1][1];
		du[2][0]    = du[2][1];
		du[0][imax] = du[0][ib2];
		du[1][imax] = du[1][ib2];
		du[2][imax] = du[2][ib2];
		
/**** left/right state with kappa = 0 ****/ 

	for(int j=1; j<=ib2; j++){
		
		minmodlim[0][j]=MinMod(du[0][j+1], du[0][j], du[0][j-1]);
		minmodlim[1][j]=MinMod(du[1][j+1], du[1][j], du[1][j-1]);
		minmodlim[2][j]=MinMod(du[2][j+1], du[2][j], du[2][j-1]);
	}
}

/****function to compute minmod limiter****/
template <typename T>
T MinMod(T Ujp2, T Ujp1, T Ujp){

	double r1,r2,del=0.000001;

	if (fabs(Ujp1)<del) r1=Sign(Ujp1)*(Ujp)/del;
	else r1=(Ujp)/(Ujp1);

	if (fabs(Ujp1)<del) r2=Sign(Ujp1)*(Ujp2)/del;
	else r2=(Ujp2)/(Ujp1);
	
	if (r1*r2<=0) return 0.0;
	else if ((fabs(r1)<fabs(r2)) and (fabs(r1)<1.0) and r1*r2>0.0) return r1;
	else if ((fabs(r2)<fabs(r1)) and (fabs(r2)<1.0) and r1*r2>0.0) return r2;
	else return 1.0;

}

template <typename T>
T Sign(T a){
	
	if (fabs(a)<1e-8) return 0;
	else if(a>0.0) return 1.0;
	else if (a<0.0) return -1.0;
	else return 0.0;
}
