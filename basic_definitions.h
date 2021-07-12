#ifndef _Header_H
#define _Header_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <string>	
#include <iosfwd>	


/***********constants used throught the program ***************/
const double cpgas=1005.0, Gamma=1.4, gam1=Gamma-1.0, gap1=Gamma+1.0, rgas=gam1*cpgas/Gamma;
const int idim=300, noconv=3, mxdum=4*idim;

using namespace std;

/***************global variables*************/
extern double **cv, **cvold, **diss, **rhs, **ls, **rs, **minmodlim, **cvrs, **cvls;
extern double *a, *p, *x, *vol, *dt, *dum;
extern int maxiter, imax,  ib2, ncells, iorder;
extern char str;
extern char *c;	
extern std::string skip, fngrid, fnplot, conver_r, rest_r;
extern double p01, t01, p2, convtol, cfl, epsirs, vis2, vis4, limfac, epsentr, nrk;
extern double *ark, *betrk, *ldiss, *lsmoo;
extern double volref, rhoref, uref, pref;



/*********************forward declaration of functions called in the main program **********************/

double **Allocate_2D(double ** &m, int t1, int t2);
void Allocate1D_Array();
void Allocate2D_Array();
void Read_input(std::ifstream & infile);
void Rgrid(ifstream & inflie,  int idim, int &imax, int &ib2, double *& x, double *& a);
void Iniflow(int &imax, double *&a, double **&cv, double *&p);
void Solver(ofstream &ofile_conv, int &iter, int &imax, int &ib2, int mxdum, double *&x, double *&a, double *&vol, double **&cv, double *&p, double **&cvold, double **&diss, double **&rhs, double *&dt, 	      double *&dum);
void Bcond(int &imax, int &ib2, double *& a, double **& cv, double *& p);
void Inigrid(int &imax, int &ib2, double *&x, double *&a, double *&vol);
void Tstep(int &imax, int &ib2, double *&x, double *&a, double *&vol, double **&cv, double *&p, double *&dt);
void LR_State_roe(int &imax, int &ib2, double *&a, double *&vol, double **&cv, double *&p);
void LR_State_movers(int &imax, int &ib2, double *&a, double *&vol, double **&cv, double *&p);
void Flux_roe(int &imax, int &ib2, double *&a, double **&ls, double **&rs);
void Flux_movers(int &imax, int &ib2, double *&a, double **&ls, double **&rs, double **cvrs, double **cvls);
void Srcterm(int &imax, int &ib2, double *&a, double *&p);
void Minmod_lim(int &imax, int &ib2, double *&a, double *&vol, double **&cv);
void Conver(ofstream &ofile, int &iter, int &imax, int &ib2, double *&a, double **&cv, double **&cvold, double *&p);
void Wsolut(int &imax, int &ib2, double *&x, double *&a, double **&cv, double *&p);




#endif  //#ifndef _Header_H
