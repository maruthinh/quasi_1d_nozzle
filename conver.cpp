#include "basic_definitions.h"

void Conver(ofstream &ofile_conv, int &iter, int &imax, int &ib2, double *&a, double **&cv, double **&cvold, double *&p){
	
	int idrho, nsup;
	double dr, drho, drho1, drmax, rho, u, c, avms;
	
/**** get the change of density and mass flow****/
	drho  = 0.0;
	drmax = -1e20;
	avms  = 0.0;
	nsup  = 0;

	for(int j=2; j<=ib2; j++){
		dr   = cv[0][j] - cvold[0][j];
		avms = avms+cv[1][j];
		drho = drho+dr*dr;
		
		if(fabs(dr)>=drmax){
			drmax=fabs(dr);
			idrho=j;
		}
		
		rho = cv[0][j]/a[j];
		u   = cv[1][j]/cv[0][j];
		c   = sqrt(Gamma*p[j]/rho);
		
		if(u>c)	nsup=nsup+1;	
	}	
	
	avms = avms/(ncells+1.0);
		
	if(iter==1){
		drho1 = sqrt(drho/(ncells+1.0));
      		drho  = sqrt(drho/(ncells+1.0))/drho1;

	}
		
	cout<<iter<<"\t"<<log10(drho)<<"\t"<<drmax<<"\t"<<idrho<<"\t"<<avms<<"\t"<<cv[1][2]-cv[1][ib2]<<"\t"<<nsup<<endl;

	ofile_conv.flags(ios::dec | ios::scientific);
	ofile_conv.precision(5);

	
	if(!ofile_conv){
		cerr<<"File couldn't be opened to write the solution"<<endl;
		exit(1);
	}

	ofile_conv<<iter<<"\t"<<drho<<"\t"<<drmax<<"\t"<<idrho<<"\t"<<avms<<"\t"<<cv[1][2]-cv[1][ib2]<<"\t"<<nsup<<endl;
	
	//ofile.close();

}
