//---This is the input.h file.
#pragma once

#ifndef INPUT_H
#define INPUT_H

//---Standard libraries
#include<cmath>
#include<vector>
#include <deque>
//---User-defined libraries


//---Namespace input

namespace input{
	
//---Variables


//---Global constants
const double mu0=4.0*M_PI*1e-7; //Vacuum permeability: [H/m]
const double gamma = 1.76*1e11; // Gyromagnetic ratio: [rad/(s*T)]
const double k_B = 1.381*1e-23; // Boltzmann constant [J/K]


//---Material parameters
extern int n_materials;
extern double lengthscale;
extern std::vector<unsigned long long int>length;
extern std::vector<unsigned long long int>width;
extern std::vector<unsigned long long int>height;
extern std::vector<unsigned long long int>unitcell_size; 
extern std::vector<double>unitcell_volume;
extern std::vector<unsigned long long int>macrocell_size; 
extern std::vector<double>macrocell_volume;
extern std::vector<unsigned long long int>material_total_cells; 
extern std::vector<unsigned long long int>material_nx_cells;
extern std::vector<unsigned long long int>material_ny_cells;
extern std::vector<unsigned long long int>material_nz_cells;
extern int nx_cells;
extern int ny_cells;
extern int nz_cells;
extern int n_cells;
extern int n_at;

extern std::vector<double> Tc;
extern std::vector<double> Ms0_CGS;
extern std::vector<double> K0_CGS;
extern std::vector<double> Ms0_SI;
extern std::vector<double> K0_SI;
extern std::vector<double> A0_list; 
extern std::vector<std::vector<double>> A0_matrix;
extern std::vector<std::vector<double>> A_T_matrix;
extern std::vector<double> m_e;
extern std::vector<double> Ms_T;
extern std::vector<double> K_T;
extern std::vector<double> mu_s;
extern std::vector<double> lambda;
extern std::vector<double> chi_par;
extern std::vector<double> chi_perp;
extern std::vector<double> eps;
extern std::vector<double> ex;
extern std::vector<double> ey;
extern std::vector<double> ez;
//---Initial conditions
extern std::vector<double>mx_0;
extern std::vector<double>my_0;
extern std::vector<double>mz_0;
extern bool read_config_from_file;

//---Simulation paramters
extern double T; // Magnetic moment temperature (electron temperature in this model): [K]
extern double B_app; // External field amplitude: [T]
extern double B_theta;//angle with respect to EA
extern double B_phi; //Angle between in-plane field projection and the Ox axis
extern double bx;// H_ext_x/H_ext 
extern double by; // H_ext_y/H_ext
extern double bz; // H_ext_z/H_ext
extern std::vector<double>alpha_par; //Longitudinal damping parameter: adimensional. Avail. only for T<Tc
extern std::vector<double>alpha_perp; //Transversal damping parameter: adimensional. Avail. only for T<Tc

extern int t_min_equil; // Time will be given as an integer. It needs to be multiplied with 10^-12 at the end.
extern int t_max_equil;
extern int delta_t_equil;
extern double timescale_equil;
extern double TOL_EQ;//Tolerance for stopping condition using the torque modulus.


extern int t_min_laser_dynamics;
extern int t_max_laser_dynamics;
extern int delta_t_laser_dynamics;
extern double timescale_laser_dynamics;
extern double TOL_LD;
extern double T_pulse;
extern int pulse_duration; //fs !!Important, this needs to be in the same time unit as the timescale_laser_dynamics

//---Simulation

extern bool m_vs_T_curve;
extern bool chipar_vs_T_curve;
extern bool K_vs_T_curve;
extern bool A_vs_T_curve;
extern bool equilibrate;
extern bool laser_dynamics;
extern bool force_DW_formation;
extern bool remove_precession_term;


//---Functions

int read_simulation_parameters();
int read_material_parameters();
int generate_exchange_matrix(int n_materials,
                              std::vector<double>A0_list,
                              std::vector<std::vector<double> >&A0_matrix);

}

#endif
//---End of input.h file.
