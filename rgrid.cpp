#include "basic_definitions.h"

/**initialization of extern global variables declared in the header file**/
double **cv, **cvold, **diss, **rhs, **ls, **rs, **minmodlim, **cvrs, **cvls;
double *a, *p, *x, *vol, *dt, *dum;
int maxiter, imax, ib2, ncells, iorder;
char str;
char *c;	
std::string skip, fngrid, fnplot, conver_r, rest_r;
double p01, t01, p2, convtol, cfl, epsirs, vis2, vis4, limfac, epsentr, nrk;
double *ark, *betrk, *ldiss, *lsmoo;
double volref, rhoref, uref, pref;

/***function to read grid co-ordinates and area from the file***/ 
void Rgrid(ifstream & infile, int idim, int &imax, int &ib2, double *& x, double *& a){

/* reads the x-coordinate and the area of the nozzle */

	//cout << "Reading co-ordintes of the nozzle and its area" << endl;
	string skip;
	infile.open("./nozzle_grid.dat");

	if(infile.fail()){
		cerr<<"File couldn't be opened to read the input parameters"<<endl;
		exit(1);
	}

	infile>>imax;

	while(!infile.eof()){
		for(int j=1;j<=imax;j++){
			infile>>x[j]>>a[j];
		}
	}

	ib2=imax-1;
	ncells=imax-3;	
}
