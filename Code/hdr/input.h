//---This is the input.h file.
#pragma once

#ifndef INPUT_H
#define INPUT_H

//---Standard libraries
#include<cmath>

//---User-defined libraries


//---Namespace input

namespace input{
	
//---Variables


//---Global constants
const double mu0=4.0*M_PI*1e-7; //Vacuum permeability: [H/m]
const double gamma = 1.76*1e11; // Gyromagnetic ratio: [rad/(s*T)]
const double k_B = 1.381*1e-23; // Boltzmann constant [J/K]


//---Material parameters
extern int n_at; //Number of atoms per unit cell: adim.
extern double a; // Lattice spacing: [m]
extern double Tc; // Curie temperature: [K]
extern double Ms0_CGS; // Saturation magnetisation: [emu/cc]
extern double K0_CGS; // Magnetocrystlline first-order anis. constant: [erg/cc]
extern double Ms0_SI; //Saturation magnetisation: [A/m]
extern double K0_SI; //: Magnetocrystlline first-order anis. constant: [J/m^3]
extern double m_e; //Equilibrium magnetisation: []
extern double Ms_T; //Saturation magnetisation at T
extern double K_T; //Magnetocrystalline anisotropy constant at T
extern double volume; // Macrospin volume: [m^3]
extern double mu_s; // Atomic magnetic moment: [J/T]
extern double lambda; // Microscopic coupling constant: adimensional
extern double chi_par; //Parallel susceptibility: []
extern double chi_perp; //Perpendicular susceptibility: []
extern double eps; //SW correction factor: adim.
extern double ex; // H_ani_x
extern double ey; // H_ani_y
extern double ez; // H_ani_z

//---Simulation paramters
extern double T; // Magnetic moment temperature (electron temperature in this model): [K]
extern double B_app; // External field amplitude: [T]
extern double B_theta;//angle with respect to EA
extern double B_phi; //Angle between in-plane field projection and the Ox axis
extern double bx;// H_ext_x/H_ext 
extern double by; // H_ext_y/H_ext
extern double bz; // H_ext_z/H_ext
extern double alpha_par; //Longitudinal damping parameter: adimensional. Avail. only for T<Tc
extern double alpha_perp; //Transversal damping parameter: adimensional. Avail. only for T<Tc

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


//---Initial conditions
extern double mx_0;
extern double my_0;
extern double mz_0;

//---Simulation

extern bool m_vs_T_curve;
extern bool chipar_vs_T_curve;
extern bool K_vs_T_curve;
extern bool equilibrate;
extern bool laser_dynamics;


//---Functions

int read_simulation_parameters();
int read_material_parameters();

}

#endif
//---End of input.h file.
