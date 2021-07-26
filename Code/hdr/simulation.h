//---This is the simulation.h file.
#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

//---Standard libraries
#include<fstream>

//---User-defined libraries


//---Namespace simulation

namespace simulation{
	
//---Variables

extern double laser_sim_time; //This holds the total simulation time of the laser-induced dynamics
extern double equil_sim_time; //This holds the total simulatino time of the equilibration process
extern double total_sim_time; //Equilibration time plus laser induced dynamics

//---Functions
int squared_pulse_dynamics(double mx_0, double my_0, double mz_0,
					   	   double &mx_n1, double &my_n1, double &mz_n1,
					   	   double &LD_REAL_t, double EQ_REAL_t,
					   	   double &T, double TOL,
					   	   int t_start, int t_end, int t_step, double timescale,
					   	   int pulse_duration, double T_pulse,
					   	   std::ofstream &f1);

int equilibrate_system(double mx_0, double my_0, double mz_0,
					   double &mx_n1, double &my_n1, double &mz_n1,
					   double &EQ_REAL_t, double T, double TOL,
					   int t_start, int t_end, int t_step, double timescale,
					   std::ofstream &f1);

}

#endif
//---End of simulation.h file.