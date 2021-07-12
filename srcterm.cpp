#include "basic_definitions.h"

void Srcterm(int &imax, int &ib2, double *&a, double *&p){
	double da;

	for(int j=2;j<=ib2;j++){
		da        = 0.5*(a[j+1]-a[j-1]);
		rhs[1][j] = rhs[1][j]-p[j]*da; 
	}
	
}
