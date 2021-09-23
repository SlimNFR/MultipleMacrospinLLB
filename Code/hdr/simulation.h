//---This is the simulation.h file.
#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

//---Standard libraries
#include<fstream>
#include<vector>
//---User-defined libraries


//---Namespace simulation

namespace simulation{
	
//---Variables

extern double laser_sim_time; //This holds the total simulation time of the laser-induced dynamics
extern double equil_sim_time; //This holds the total simulatino time of the equilibration process
extern double total_sim_time; //Equilibration time plus laser induced dynamics

//---Functions
int force_DW_formation_f(std::vector<double> m_e, std::vector<double> &mx, std::vector<double> &my, std::vector<double> &mz,
						 std::vector<int>material_id);

int squared_pulse_dynamics(int n_cells,
						   double gamma, std::vector<double> &alpha_par, std::vector<double>&alpha_perp,
						   std::vector<double> mx_0, std::vector<double> my_0, std::vector<double> mz_0,
						   std::vector<double> &mx_n1, std::vector<double> &my_n1, std::vector<double> &mz_n1,
						   std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff,
						   std::vector<double> &torque_mod,
						   std::vector<int> material_id,
					   	   double &LD_REAL_t, double EQ_REAL_t,
					   	   double &T, double TOL,
					   	   int t_start, int t_end, int t_step, double timescale,
					   	   int pulse_duration, double T_pulse,
					   	   std::ofstream &f1,
					   	   std::ofstream &f2);

int equilibrate_system(int n_cells,
					   double gamma, std::vector<double> &alpha_par, std::vector<double>&alpha_perp,
					   std::vector<double> mx_0, std::vector<double> my_0, std::vector<double> mz_0,
					   std::vector<double> &mx_n1, std::vector<double> &my_n1, std::vector<double> &mz_n1,
					   std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff,
					   std::vector<double> &torque_mod,
					   std::vector<int> material_id,
					   double &EQ_REAL_t, double T, double TOL,
					   int t_start, int t_end, int t_step, double timescale,
					   std::ofstream &f1,
					   std::ofstream &f2);

}

#endif
//---End of simulation.h file.