#include "basic_definitions.h"

/**** Forward declaration of functions****/
template <typename T>
T MinModLim(T Ur, T Ul);


/****computes left and right state using MUSCL appraoch****/

void LR_State_movers(int &imax, int &ib2, double *&a, double *&vol, double **&cv, double *&p){
	
	double **du;
	double *dx, *deltl, *deltr;
	
/****dynamic memory allocation for 2D & 1D pointers****/
	Allocate_2D(du, noconv, idim);
	dx=new double[idim]; deltl=new double[3]; deltr=new double[3];
		
/****if iorder==1, then first order, iorder==2, for second ****/
			
	

	if (iorder==2){
/****first difference of rho, u, p****/
		for(int j=1;j<=ib2;j++){
			du[0][j]    = cv[0][j+1] - cv[0][j];
			du[1][j]    = cv[1][j+1] - cv[1][j];
			du[2][j]    = cv[2][j+1] - cv[2][j];
		}
	
		du[0][0]    = du[0][1];
		du[1][0]    = du[1][1];
		du[2][0]    = du[2][1];
		du[0][imax] = du[0][ib2];
		du[1][imax] = du[1][ib2];
		du[2][imax] = du[2][ib2];
		
		for(int j=2;j<=ib2;j++){
			dx[j] = 0.5*(x[j+1]-x[j-1]);
			dx[j]=sqrt(dx[j]*dx[j]);
		}
		
		dx[1] = dx[2];
		dx[imax]=dx[ib2];

 		for(int j=1; j<=ib2; j++){
			deltr[0] = 0.5*MinModLim(du[0][j+1], du[0][j]);
			deltl[0] = 0.5*MinModLim(du[0][j],   du[0][j-1]);
			deltr[1] = 0.5*MinModLim(du[1][j+1], du[1][j]);
			deltl[1] = 0.5*MinModLim(du[1][j],   du[1][j-1]);
			deltr[2] = 0.5*MinModLim(du[2][j+1], du[2][j]);
			deltl[2] = 0.5*MinModLim(du[2][j],   du[2][j-1]);
			cvrs[0][j] = cv[0][j+1] - deltr[0];
			cvrs[1][j] = cv[1][j+1] - deltr[1];
			cvrs[2][j] = cv[2][j+1] - deltr[2];
			cvls[0][j] = cv[0][j]   + deltl[0];
			cvls[1][j] = cv[1][j]   + deltl[1];
			cvls[2][j] = cv[2][j]   + deltl[2];
		}

		for(int j=1; j<=ib2; j++){
			rs[0][j] = cv[0][j+1]/a[j+1];
			rs[1][j] = cv[1][j+1]/cv[0][j+1];
			rs[2][j] = p[j+1];
			ls[0][j] = cv[0][j]/a[j];
			ls[1][j] = cv[1][j]/cv[0][j];
			ls[2][j] = p[j];
		}
	}

/**** for first order ****/
	else if(iorder==1){
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
		cout<<j<<"\t"<<rs[1][j]<<"\t"<<rs[2][j]<<"\t"<<rs[3][j]<<"\t"<<ls[1][j]<<"\t"<<ls[2][j]<<"\t"<<ls[3][j]<<endl;
	}*/
	
	//cout<<"reference values"<<limfac3<<"\t"<<rvolref<<"\t"<<eps2[0]<<"\t"<<eps2[1]<<"\t"<<eps2[2]<<endl;
	//cout<<"akkan iorder"<<eps2[0]<<"\t"<<eps2[1]<<"\t"<<eps2[2]<<"\t"<<vola<<"\t"<<eps2n<<"\t"<<deltr[1]<<endl;
}	


/****Minmod slope limiter***/
template <typename T>
T MinModLim(T Ur, T Ul){
	
	double Phi;
	if ((fabs(Ur)<fabs(Ul)) and Ur*Ul>0) Phi = Ur;
	else if ((fabs(Ul)<fabs(Ur)) and Ur*Ul>0) Phi = Ul;
	else if (Ur*Ul<=0) Phi = 0.0;
	
	return Phi;
}

