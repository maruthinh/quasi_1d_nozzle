#include "basic_definitions.h"


/****function to compute control volumes****/

void Inigrid(int &imax, int &ib2, double *&x, double *&a, double *&vol){
	
	double da;

	for(int j=2; j<=ib2; j++){
		da=a[j-1]+2.0*a[j]+a[j+1];
		vol[j]=(x[j+1]-x[j-1])*da/8.0;
	}
	
	vol[1]=vol[2];
	vol[imax]=vol[ib2];

}
