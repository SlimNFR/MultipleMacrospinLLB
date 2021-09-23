//---This is the solver.cpp file. It defines the solver.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>

//---User-defined libraries
#include"solver.h"
//---Namespace solver

namespace solver{
	
//---Variables

int sim_time = 0;

//---Functions

int heun_scheme_step(int (*dfunc)(bool,
					 double, double, double,
					 double, double, double,
					 double, double, double,
					 double &,double &,double &),
					 bool remove_precession_term,
					 double mx_0,double my_0, double mz_0,
					 double gamma, double alpha_par, double alpha_perp,
					 int delta_t, double timescale, //timescale variable should be intialised to 1ps or 1ns or 1fs depending on the timescale of my sim
					 double Bx_eff, double By_eff, double Bz_eff,
					 double &mx_n1, double &my_n1, double &mz_n1)
{	//This function returns the next magnetisatino step given initial conditions and parameters for the LLB equation.
	double fx_0,fy_0,fz_0 = 0.0;
	double fx_n,fy_n,fz_n = 0.0;
	double REAL_DT = delta_t*timescale;// the for loop imaginary time step transformed in ps;
		/*
		std::cout<<"|fx_0: "<<fx_0<<"|fy_0: "<<fy_0<<"|fz_0: "<<fz_0
				 <<"|fx_n: "<<fx_n<<"|fy_n: "<<fy_n<<"|fz_n: "<<fz_n<<"\n";*/
		/*std::cout<<REAL_T<<": |"<<"Bx_eff: "<<field::Bx_eff<<"|By_eff: "<<field::By_eff<<"|Bz_eff: "<<field::Bz_eff
		<<"|Bx_ani: "<<field::Bx_ani<<"|By_ani: "<<field::By_ani<<"|Bz_ani: "<<field::Bz_ani
		<<"|mx: "<<mx_n1<<"|my: "<<my_n1<<"|mz: "<<mz_n1<<"\n";
*/		
		//std::cout<<"|gamma: "<<gamma<<"|alpha_par: "<<alpha_par<<"|alpha_perp: "<<alpha_perp<<"\n";
		dfunc(remove_precession_term,mx_0,my_0,mz_0,Bx_eff,By_eff,Bz_eff,gamma,alpha_par,alpha_perp,fx_0,fy_0,fz_0);
		std::cout<<fx_0<<" "<<fy_0<<" "<<fz_0<<"\n";
		dfunc(remove_precession_term,mx_0+REAL_DT*fx_0,my_0+REAL_DT*fy_0,mz_0+REAL_DT*fz_0,Bx_eff,By_eff,Bz_eff,gamma,alpha_par,alpha_perp,fx_n,fy_n,fz_n);

		mx_n1 = mx_0  + (REAL_DT/2.0)*(fx_0 + fx_n);
		my_n1 = my_0  + (REAL_DT/2.0)*(fy_0 + fy_n);
		mz_n1 = mz_0  + (REAL_DT/2.0)*(fz_0 + fz_n);

		/*std::cout<<"|mx: "<<mx_n1<<"|my: "<<my_n1<<"|mz: "<<mz_n1
				 <<"|mx_0: "<<mx_0<<"|my_0: "<<my_0<<"|mz_0: "<<mz_0
				 <<"|fx_0 + fx_n: "<<fx_0 + fx_n<<"|fy_0 + fy_n: "<<fy_0 + fy_n<<"|fz_0 + fz_n: "<<fz_0 + fz_n
				 <<"|REAL_DT: "<<REAL_DT<<"\n";
		*/

	return 0;

}


int RK4_scheme_step(int (*dfunc)(bool,
					 double, double, double,
					 double, double, double,
					 double, double, double,
					 double &,double &,double &),
					 bool remove_precession_term,
					 double mx_0,double my_0, double mz_0,
					 double gamma, double alpha_par, double alpha_perp,
					 int delta_t, double timescale, //timescale variable should be intialised to 1ps or 1ns or 1fs depending on the timescale of my sim
					 double Bx_eff, double By_eff, double Bz_eff,
					 double &mx_n1, double &my_n1, double &mz_n1)
{	//This function returns the next magnetisatino step given initial conditions and parameters for the LLB equation.
	double REAL_DT = delta_t*timescale;// the for loop imaginary time step transformed in ps;


	double k1_x, k1_y, k1_z = 0.0;
	double k2_x, k2_y, k2_z = 0.0;
	double k3_x, k3_y, k3_z = 0.0;
	double k4_x, k4_y, k4_z = 0.0;
	double k_x, k_y, k_z = 0.0;



	//Calculate k1
	dfunc(remove_precession_term,mx_0,my_0,mz_0,Bx_eff,By_eff,Bz_eff,gamma,alpha_par,alpha_perp,k1_x,k1_y,k1_z);
	//std::cout<<k1_x<<" "<<k1_y<<" "<<k1_z<<"\n";		
	//Calculate k2
	dfunc(remove_precession_term,mx_0+(REAL_DT/2.0)*k1_x,my_0+(REAL_DT/2.0)*k1_y,mz_0+(REAL_DT/2.0)*k1_z,Bx_eff,By_eff,Bz_eff,gamma,alpha_par,alpha_perp,k2_x,k2_y,k2_z);
	//std::cout<<k2_x<<" "<<k2_y<<" "<<k2_z<<"\n";	
	//Calculate k3
	dfunc(remove_precession_term,mx_0+(REAL_DT/2.0)*k2_x,my_0+(REAL_DT/2.0)*k2_y,mz_0+(REAL_DT/2.0)*k2_z,Bx_eff,By_eff,Bz_eff,gamma,alpha_par,alpha_perp,k3_x,k3_y,k3_z);
	//std::cout<<k3_x<<" "<<k3_y<<" "<<k3_z<<"\n";	
	//Calculate k4
	dfunc(remove_precession_term,mx_0+REAL_DT*k3_x,my_0+REAL_DT*k3_y,mz_0+REAL_DT*k3_z,Bx_eff,By_eff,Bz_eff,gamma,alpha_par,alpha_perp,k4_x,k4_y,k4_z);	
	//std::cout<<k4_x<<" "<<k4_y<<" "<<k4_z<<"\n";	
	//Calculate k


	k_x = (1.0/6.0)*REAL_DT*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
	k_y = (1.0/6.0)*REAL_DT*(k1_y + 2.0*k2_y + 2.0*k3_y + k4_y);
	k_z = (1.0/6.0)*REAL_DT*(k1_z + 2.0*k2_z + 2.0*k3_z + k4_z);

		mx_n1 = mx_0  + k_x;
		my_n1 = my_0  + k_y;
		mz_n1 = mz_0  + k_z;


	return 0;

}


}


//---End of solver.cpp file.

