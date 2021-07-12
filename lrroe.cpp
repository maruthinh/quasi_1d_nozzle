#include "basic_definitions.h"

/**** Forward declaration of functions****/
double MUSCL0(double af, double bf, double epsilon);
double MUSCL3(double af, double bf, double epsilon);


/****computes left and right state using MUSCL appraoch****/

void LR_State_roe(int &imax, int &ib2, double *&a, double *&vol, double **&cv, double *&p){
	
	double **du;
	double *eps2, *deltl, *deltr;
	double limfac3, rvolref, vola, eps2n;

/****dynamic memory allocation for 2D & 1D pointers****/
	Allocate_2D(du, noconv, idim);
	eps2=new double[3]; deltl=new double[3]; deltr=new double[3];
		
/****if iorder==1, then first order, iorder==2, for kappa=0, and iorder==3, for kappa=1/3 ****/
			
	
	if (iorder==2 or iorder==3){
/****normalised epsilon^2 for all limited variables (rho, u, p)****/
		limfac3 = limfac*limfac*limfac;
		rvolref = 1.0/pow(volref,1.5);
		eps2[0]  = limfac3*rhoref*rhoref*rvolref;
		eps2[1]  = limfac3*uref*uref*rvolref;
		eps2[2]  = limfac3*pref*pref*rvolref;

/****first difference of rho, u, p****/
		for(int j=1;j<=ib2;j++){
			du[0][j]    = cv[0][j+1]/a[j+1] - cv[0][j]/a[j];
			du[1][j]    = cv[1][j+1]/cv[0][j+1] - cv[1][j]/cv[0][j];
			du[2][j]    = p[j+1]-p[j];
		}
	
		du[0][0]    = du[0][1];
		du[1][0]    = du[1][1];
		du[2][0]    = du[2][1];
		du[0][imax] = du[0][ib2];
		du[1][imax] = du[1][ib2];
		du[2][imax] = du[2][ib2];

		if (iorder==2 and iorder!=3){

/**** left/right state with kappa = 0 ****/ 
			for(int j=1; j<=ib2; j++){
				vola     = pow(0.5*(vol[j+1]+vol[j]), 1.5);
				eps2n    = eps2[0]*vola;
				deltr[0] = 0.5*MUSCL0(du[0][j+1], du[0][j], eps2n);
				deltl[0] = 0.5*MUSCL0(du[0][j], du[0][j-1], eps2n);
				eps2n    = eps2[1]*vola;
				deltr[1] = 0.5*MUSCL0(du[1][j+1], du[1][j], eps2n);
				deltl[1] = 0.5*MUSCL0(du[1][j], du[1][j-1], eps2n);
				eps2n    = eps2[2]*vola;
				deltr[2] = 0.5*MUSCL0(du[2][j+1], du[2][j], eps2n);
				deltl[2] = 0.5*MUSCL0(du[2][j], du[2][j-1], eps2n);
				rs[0][j] = cv[0][j+1]/a[j+1] 	 - deltr[0];
				rs[1][j] = cv[1][j+1]/cv[0][j+1] - deltr[1];
				rs[2][j] = p[j+1] 	         - deltr[2];
				ls[0][j] = cv[0][j]/a[j] 	 + deltl[0];
				ls[1][j] = cv[1][j]/cv[0][j]     + deltl[1];
				ls[2][j] = p[j]        		 + deltl[2];
			}
		}

		else if(iorder==3 and iorder!=2){
/**** left/right state with kappa = 1/3 ****/ 
			for(int j=1; j<=ib2; j++){
				vola     = pow((0.5*(vol[j+1]+vol[j])), 1.5);
				eps2n    = eps2[0]*vola;
				deltr[0] = 0.5*MUSCL3(du[0][j+1], du[0][j], eps2n);
				deltl[0] = 0.5*MUSCL3(du[0][j], du[0][j-1], eps2n);
				eps2n    = eps2[1]*vola;
				deltr[1] = 0.5*MUSCL3(du[1][j+1], du[1][j], eps2n);
				deltl[1] = 0.5*MUSCL3(du[1][j], du[1][j-1], eps2n);
				eps2n    = eps2[2]*vola;
				deltr[2] = 0.5*MUSCL3(du[2][j+1], du[2][j], eps2n);
				deltl[2] = 0.5*MUSCL3(du[2][j], du[2][j-1], eps2n);
				rs[0][j] = cv[0][j+1]/a[j+1] 	 - deltr[0];
				rs[1][j] = cv[1][j+1]/cv[0][j+1] - deltr[1];
				rs[2][j] = p[j+1] 	         - deltr[2];
				ls[0][j] = cv[0][j]/a[j] 	 + deltl[0];
				ls[1][j] = cv[1][j]/cv[0][j]     + deltl[1];
				ls[2][j] = p[j]        		 + deltl[2];
			}	
		}
	}
/**** for first order ****/
	else if(iorder==1 and iorder!=2 and iorder!=3){
		for(int j=1; j<=ib2; j++){
			rs[0][j] = cv[0][j+1]/a[j+1];
			rs[1][j] = cv[1][j+1]/cv[0][j+1];
			rs[2][j] = p[j+1];
			ls[0][j] = cv[0][j]/a[j];
			ls[1][j] = cv[1][j]/cv[0][j];
			ls[2][j] = p[j];
		}

	}
	else{
		cout<<"wrong input for order of accuracy"<<endl;
		exit(0);
	}			
	/*for(int j=1; j<=imax; j++){
		cout<<j<<"\t"<<rs[0][j]<<"\t"<<rs[1][j]<<"\t"<<rs[2][j]<<endl;
	}*/
	
	//cout<<"reference values"<<limfac3<<"\t"<<rvolref<<"\t"<<eps2[0]<<"\t"<<eps2[1]<<"\t"<<eps2[2]<<endl;
	//cout<<"akkan iorder"<<eps2[0]<<"\t"<<eps2[1]<<"\t"<<eps2[2]<<"\t"<<vola<<"\t"<<eps2n<<"\t"<<deltr[1]<<endl;
}	


/****MUSCL scheme with kappa=0***/
double MUSCL0(double af, double bf, double epsilon){

	return (af*(bf*bf+epsilon)+bf*(af*af+epsilon))/(af*af+bf*bf+2.0*epsilon+1e-10);
}

/****MUSCL scheme with kappa=1/3***/
double MUSCL3(double af, double bf, double epsilon){

	return (bf*(2.0*af*af+epsilon)+af*(bf*bf+2.0*epsilon))/(2.0*af*af+2.0*bf*bf-af*bf+3.0*epsilon+1e-10);
}

