//---This is the solver.cpp file. It defines the solver.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>
#include<vector>

//---User-defined libraries
#include"solver.h"
#include"structure.h"
#include"field.h"
//---Namespace solver

namespace solver{
	
//---Variables

int sim_time = 0;

//---Functions

int heun_scheme_step(int n_cells,
					 int (*dfunc)(bool, bool, bool,
					 double, double, double,
					 double, double, double,
					 double, double, double,
					 double &,double &,double &),
					 bool remove_precession_term,
					 bool remove_longitudin_term,
					 bool remove_transverse_term,
					 std::vector<int>material_id,
					 std::vector<double> &mx_0, std::vector<double> &my_0, std::vector<double> &mz_0,
					 std::vector<double> &mx_n1, std::vector<double> &my_n1, std::vector<double> &mz_n1,
					 std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff,
					 double gamma, std::vector<double> alpha_par, std::vector<double> alpha_perp,
					 int delta_t, double timescale
					 )//timescale variable should be intialised to 1ps or 1ns or 1fs depending on the timescale of my sim
					 
{	//This function returns the next magnetisatino step given initial conditions and parameters for the LLB equation.
	std::vector<double> fx_0, fy_0, fz_0;
	std::vector<double> fx_n, fy_n, fz_n;

	fx_0.resize(n_cells,0.0);
	fy_0.resize(n_cells,0.0);
	fz_0.resize(n_cells,0.0);

	fx_n.resize(n_cells,0.0);
	fy_n.resize(n_cells,0.0);
	fz_n.resize(n_cells,0.0);


	double REAL_DT = delta_t*timescale;


	//First Heun step
	for(int cell=0; cell<n_cells; cell++)
	{
	//
		int mat_id = material_id[cell];
		dfunc(remove_precession_term,remove_longitudin_term,remove_transverse_term,
			  mx_0[cell],my_0[cell],mz_0[cell],
			  Bx_eff[cell],By_eff[cell],Bz_eff[cell],
			  gamma,alpha_par[mat_id],alpha_perp[mat_id],
			  fx_0[cell],fy_0[cell],fz_0[cell]);
	//
		macrospin::mx[cell]=mx_0[cell]+REAL_DT*fx_0[cell];
		macrospin::my[cell]=my_0[cell]+REAL_DT*fy_0[cell];
		macrospin::mz[cell]=mz_0[cell]+REAL_DT*fz_0[cell];

	}
	

	field::calculate();

	//Second Heun step

	for(int cell=0; cell<n_cells; cell++)
	{
		int mat_id = material_id[cell];
		dfunc(remove_precession_term,remove_longitudin_term,remove_transverse_term,
			  mx_0[cell]+REAL_DT*fx_0[cell],my_0[cell]+REAL_DT*fy_0[cell],mz_0[cell]+REAL_DT*fz_0[cell],
			  Bx_eff[cell],By_eff[cell],Bz_eff[cell],
			  gamma,alpha_par[mat_id],alpha_perp[mat_id],
			  fx_n[cell],fy_n[cell],fz_n[cell]);
		
	}

	//Update macrospin vectors and set new initial conditions.
	for(int cell=0; cell<n_cells; cell++)
	{
		int mat_id = material_id[cell];
		mx_n1[cell] = mx_0[cell]  + (REAL_DT/2.0)*(fx_0[cell] + fx_n[cell]);  
		my_n1[cell] = my_0[cell]  + (REAL_DT/2.0)*(fy_0[cell] + fy_n[cell]);
		mz_n1[cell] = mz_0[cell]  + (REAL_DT/2.0)*(fz_0[cell] + fz_n[cell]);



		mx_0[cell] = mx_n1[cell];
		my_0[cell] = my_n1[cell];
		mz_0[cell] = mz_n1[cell];
	}
	
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


	/*
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

	*/
	return 0;

}


}


//---End of solver.cpp file.

