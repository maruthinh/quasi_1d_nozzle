#include "basic_definitions.h"

/**************  to allocate 1d array***********/


void Allocate1D_Array(){
	a = new double[idim];  p = new double[idim];  x = new double[idim];  vol = new double[idim];
	dt = new double[idim];  dum = new double[idim];	c=new char[10]; ark = new double[10]; 
}


/**************  to allocate 2d array***********/
void Allocate2D_Array(){
	Allocate_2D(cv, noconv, idim);	 Allocate_2D(cvold, noconv, idim); Allocate_2D(diss, noconv, idim);	
	Allocate_2D(rhs, noconv, idim);	 Allocate_2D(ls, noconv, idim);	   Allocate_2D(rs, noconv, idim);	
	Allocate_2D(minmodlim, noconv, idim); Allocate_2D(cvrs, noconv, idim); Allocate_2D(cvls, noconv, idim);
}

/******this function reads the numers in line in to array by skipping the commas seperating them******/
double *Read_line(ifstream & file, int size){
	char c;
	double *array;
	array=new double[size];	

	for(int i=0;i<size+1;i++){
		file>>array[i]>>c;
	}
	file.ignore(80,'\n');
	return array;
}

/******this function reads the necessary input from the file*********/
void Read_input(ifstream & infile){

	infile.open("./input_r");
	
	if(infile.fail()){
		cerr<<"File couldn't be opened to read the input parameters"<<endl;
		exit(1);
	}

/************to read input and output files****************/

	getline(infile, skip);
	getline(infile, fngrid);
	getline(infile, skip);
	getline(infile, skip);
	getline(infile, fnplot);
	getline(infile, skip);
	getline(infile, skip);
	getline(infile, conver_r);
	getline(infile, skip);
	getline(infile, skip);
	getline(infile, rest_r);
	getline(infile, skip);
	getline(infile, skip);
	getline(infile, skip);

/*****************to read input physical quantities************/	
	infile>>p01;
	infile.ignore(80,'\n');
	infile>>t01;
	infile.ignore(80,'\n');
	infile>>p2;
	infile.ignore(80,'\n');
	getline(infile, skip);
	getline(infile, skip);
	getline(infile, skip);

/*****************to read input characters to control the iterations************/	
	infile>>maxiter;
	infile.ignore(80,'\n');
	infile>>convtol;
	infile.ignore(80,'\n');
	infile>>c[0];
	infile.ignore(80,'\n');
	getline(infile, skip);
	getline(infile, skip);
	getline(infile, skip);

/*****************to read numerical parameters************/	
	infile>>cfl;
	infile.ignore(80,'\n');
	infile>>epsirs;
	infile.ignore(80,'\n');
	infile>>c[1];
	infile.ignore(80,'\n');
	infile>>vis2;
	infile.ignore(80,'\n');
	infile>>vis4;
	infile.ignore(80,'\n');
	infile>>iorder;
	infile.ignore(80,'\n');
	infile>>limfac;
	infile.ignore(80,'\n');
	infile>>epsentr;
	infile.ignore(80,'\n');
	infile>>nrk;
	infile.ignore(80,'\n');
	ark=Read_line(infile, 4);
	betrk=Read_line(infile, 4);
	ldiss=Read_line(infile, 4);
	lsmoo=Read_line(infile, 4);
	
	

	infile.close();
	cout<<"-------------------------------------------------------------"<<endl   
	    <<" Reading input and output files to read and write " <<endl
	    <<"-------------------------------------------------------------"<<endl
	    <<" name of the grid file to read is                	    = "<<fngrid<<endl
	    <<" name of the file to write solution is           	    = "<<fnplot<<endl
	    <<" name of the file to write convergence history is	    = "<<conver_r<<endl
	    <<" name of the file to from which the solution is restarted is = "<<rest_r<<endl
	    <<"-------------------------------------------------------------"<<endl
	    <<" reading physical quantities such as total pressure and temperature..etc"<<endl
	    <<"-------------------------------------------------------------"<<endl
	    <<" total inlet pressure is    = " << p01 <<endl				
	    <<" total inlet temperature is = " << t01 <<endl
	    <<" outlet static pressure     = " << p2 <<endl
            <<"-------------------------------------------------------------"<<endl
	    <<" reading the parameters to control the iterations...."<<	endl
	    <<"-------------------------------------------------------------"<<endl
	    <<" The maximum number of iterations    = " << maxiter<<endl
	    <<" Convergence tolerence               = " << convtol <<endl
            <<" Restart the solution from the file  = " << c[0] <<endl		
	    <<"-------------------------------------------------------------"<<endl			
	    <<" Numerical parameters required for the program "<<endl	
	    <<"-------------------------------------------------------------"<<endl	
	    <<" cfl number for time-step calculation            = " << cfl <<endl 	
	    <<" coefficient of implicit residual smoothing      = " << epsirs <<endl
	    <<" Type of the scheme used                         = " << c[1] <<endl
	    <<" 2nd order artificial dissipation coefficient-k2 = " << vis2 << endl
	    <<" 4nd order artificial dissipation coefficient-k4 = " << vis4 << endl
	    <<" Spacial order for Roe's scheme                  = " << iorder << endl
	    <<" limiter coefficient (Roe scheme)                = " << limfac << endl
	    <<" Entropy correction coefficient for Roe scheme   = " << epsentr <<endl
	    <<" Number of stages for R-K method                 = " << nrk << endl
	    <<" The stage coefficient for R-K method are below, " << endl ;
	for(int i=0;i<5;i++){
	    	    cout<<ark[i]<<"\t";
	}

	cout<<endl	
	    <<" dissipation blending coefficients are  " <<endl;

	for(int i=0;i<5;i++){
	    	    cout<<betrk[i]<<"\t";
        }	

	cout<<endl	
	    <<" dissipation evaluation, 1 for yes and 0 for no  " <<endl;

	for(int i=0;i<5;i++){
	    	    cout<<ldiss[i]<<"\t";
	}

	cout<<endl	
	    <<" residual smoothing, 1 for yes and 0 for no  " <<endl;

	for(int i=0;i<5;i++){
	    	    cout<<lsmoo[i]<<"\t";
        }

	cout<<endl;
}

int main(){
	
	int iter=0;
	ifstream filename;
	ifstream file_grid;
	ofstream ofile_conv("convergence.dat");
	ofile_conv<<"step"<<"\t"<<"log(res)"<<"\t"<<"drho_max"<<"\t"<<"node"<<"\t"<<"mass_flow"<<"\t"<<"dmass"<<"\t"<<"\t"<<"nsup"<<endl;
	Allocate1D_Array();
	Allocate2D_Array();	
 	Read_input(filename);
	Rgrid(file_grid, idim, imax, ib2, x, a);
	Inigrid(imax, ib2, x, a, vol);
	Iniflow(imax, a, cv, p);
	Bcond(imax, ib2, a, cv, p);

	while (iter<=maxiter){
				
		Solver(ofile_conv, iter, imax, ib2, mxdum, x, a, vol, cv, p, cvold, diss, rhs, dt, dum);		

		iter=iter+1;
	}
	
	Wsolut(imax, ib2, x, a, cv, p);


	return 0;
		
}


