#include "basic_definitions.h"

/****function to write sollution output to file****/

void Wsolut(int &imax, int &ib2, double *&x, double *&a, double **&cv, double *&p){
	
	double rho, u, temp, c, mach;
	
	ofstream ofile("solution.dat");
	ofile.flags(ios::dec | ios::scientific);
	ofile.precision(5);

	
	if(!ofile){
		cerr<<"File couldn't be opened to write the solution"<<endl;
		exit(1);
	}

	ofile<<"x"<<"\t"<<"\t"<<"a"<<"\t"<<"\t"<<"rho"<<"\t"<<"\t"<<"u"<<"\t"<<"\t"<<"p"<<"\t"<<"\t"<<"tempr"<<"\t"<<"\t"<<"mach"<<"\t"<<"\t"<<"massflow"<<endl;

	for(int j=2; j<=ib2; j++){
		rho  = cv[0][j]/a[j];
		u    = cv[1][j]/cv[0][j];
		temp = p[j]/(rgas*rho);
		c    = sqrt(Gamma*p[j]/rho);
		mach = u/c;
	
		ofile<<x[j]<<"\t"<<a[j]<<"\t"<<rho<<"\t"<<u<<"\t"<<p[j]<<"\t"<<temp<<"\t"<<mach<<"\t"<<cv[1][j]<<endl;
	}
		
	ofile.close();
}
