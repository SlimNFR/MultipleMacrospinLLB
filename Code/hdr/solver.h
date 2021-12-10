//---This is the solver.h file.
#pragma once

#ifndef SOLVER_H
#define SOLVER_H

//---Standard libraries
#include<vector>

//---User-defined libraries


//---Namespace solver

namespace solver{
	
//---Variables
extern int sim_time;

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
					 int delta_t, double timescale);

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
					 double &mx_n1, double &my_n1, double &mz_n1);

}

#endif
//---End of solver.h file.