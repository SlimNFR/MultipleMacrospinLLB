//---This is the input.cpp file. It defines the input.h functions.


//---Standard libraries
#include<cmath>
#include<iostream>
#include<fstream>
#include<string>

//---User-defined libraries
#include"input.h"

//---Namespace input

namespace input{
	

//---Variables

// //---Global constants


// //---Material parameters
int n_at; //Number of atoms per unit cell: adim.
double a; // Lattice spacing: [m]
double Tc; // Curie temperature: [K]
double Ms0_CGS; // Saturation magnetisation: [emu/cc]
double K0_CGS; // Magnetocrystlline first-order anis. constant: [erg/cc]
double Ms0_SI; //Saturation magnetisation: [A/m]
double K0_SI; //: Magnetocrystlline first-order anis. constant: [J/m^3]
double m_e; //Equilibrium magnetisation: []
double Ms_T; //Saturation magnetisation at T
double K_T;// Magnetocrystalline anisotropy constant at T
double volume; // Macrospin volume: [m^3]
double mu_s; // Atomic magnetic moment: [J/T]
double lambda; // Microscopic coupling constant: adimensional
double chi_par; //Parallel susceptibility: []
double chi_perp; //Perpendicular susceptibility: []
double eps; //SW correction factor: adim.
double ex; // H_ani_x
double ey; // H_ani_y
double ez; // H_ani_z

// //---Simulation paramters
double T; // Magnetic moment temperature (electron temperature in this model): [K]
double B_app; // External field amplitude: [T]
double B_theta;//angle with respect to EA
double B_phi; //Angle between in-plane field projection and the Ox axis
double bx; // H_ext_x/H_ext 
double by; // H_ext_y/H_ext
double bz; // H_ext_z/H_ext
double alpha_par; //Longitudinal damping parameter: adimensional. Avail. only for T<Tc
double alpha_perp; //Transversal damping parameter: adimensional. Avail. only for T<Tc

int t_min_equil; // Time will be given as an integer. It needs to be multiplied with timescale_equil to obtain the real-time.
int t_max_equil;
int delta_t_equil;
double timescale_equil;
double TOL_EQ;//Tolerance for stopping condition using the torque modulus.


int t_min_laser_dynamics;
int t_max_laser_dynamics;
int delta_t_laser_dynamics;
double timescale_laser_dynamics;
double TOL_LD;
double T_pulse;
int pulse_duration; //fs !!Important, this needs to be in the same time unit as the timescale_laser_dynamics variable

//---Initial conditions
double mx_0;
double my_0;
double mz_0;

//---Simulation

bool m_vs_T_curve;
bool chipar_vs_T_curve;
bool K_vs_T_curve;
bool equilibrate;
bool laser_dynamics;


//---Functions
int read_simulation_parameters()
{

 std::string s1;
 std::string input_file="simInput.txt";
 std::ifstream inFile(input_file);


if(!inFile)
    {
        std::cout << std::endl << "Failed to open the file: "<< input_file <<std::endl;
        exit(1);
    }

    std::cout<<"Simulation parameters: "<<"\n";
    std::cout<<std::endl;
    inFile >> s1>> input::T;
    std::cout <<"Initial Temperature [K]:"<< " " << input::T << "\t"<< std::endl;
  	inFile >> s1>> input::B_app;
    std::cout <<"Applied field strength [T]:"<< " " << input::B_app << "\t"<< std::endl;
    inFile >> s1>> input::B_theta;
    input::B_theta=input::B_theta*M_PI/180.0;
    std::cout <<"Applied field theta angle [degr.]:"<< " " << input::B_theta << "\t"<< std::endl;
    inFile >> s1>> input::B_phi;
    input::B_phi=input::B_phi*M_PI/180.0;
    std::cout <<"Applied field phi angle [degr.]:"<< " " << input::B_phi << "\t"<< std::endl;  
    input::bx= sin(input::B_theta)*cos(input::B_phi); // H_ext_x/H_ext 
	input::by = sin(input::B_theta)*sin(input::B_phi); // H_ext_y/H_ext
	input::bz = cos(input::B_theta); // H_ext_z/H_ext
	std::cout<<"Applied field X coordinate [adim.]:"<<" "<<input::bx<<"\n";
	std::cout<<"Applied field Y coordinate [adim.]:"<<" "<<input::by<<"\n";
	std::cout<<"Applied field Z coordinate [adim.]:"<<" "<<input::bz<<"\n";
  	inFile >> s1>> input::t_min_equil;
    std::cout <<"Initial equilibration time start [adim.]:"<< " " << input::t_min_equil << "\t "<< std::endl;
    inFile >> s1>> input::t_max_equil;
    std::cout <<"Initial equilibration time stop [adim.]:"<< " " << input::t_max_equil << "\t "<< std::endl;
    inFile >> s1>> input::delta_t_equil;
    std::cout <<"Equilibration timestep [adim.]:" << " " << input::delta_t_equil << "\t"<< std::endl;
    inFile >> s1>> input::timescale_equil;
    std::cout <<"Equilibration timescale [s]:"<< " " << input::timescale_equil << "\t"<< std::endl;
    inFile >> s1>> input::TOL_EQ;
    std::cout <<"Equilibration convergence tolerance [adim.]:"<< " " << input::TOL_EQ << "\t"<< std::endl;

    inFile >> s1>> input::t_min_laser_dynamics;
    std::cout <<"Laser dynamics time start [adim.]:"<< " " << input::t_min_laser_dynamics << "\t"<< std::endl;
    inFile >> s1>> input::t_max_laser_dynamics;
    std::cout <<"Laser dynamics time start [adim.]:"<< " " << input::t_max_laser_dynamics << "\t"<< std::endl;
    inFile >> s1>> input::delta_t_laser_dynamics;
    std::cout <<"Laser dynamics timestep [adim.]:"<<" "<< input::delta_t_laser_dynamics << " "<< std::endl;
    inFile >> s1>> input::timescale_laser_dynamics;
    std::cout <<"Laser dynamics timescale [s]:"  << " " << input::timescale_laser_dynamics << "\t"<< std::endl;
    inFile >> s1>> input::TOL_LD;
    std::cout <<"Laser dynamics convergence tolerance [adim.]:"<< " " << input::TOL_LD << "\t"<< std::endl;
    inFile >> s1>> input::T_pulse;
    std::cout <<"Temperature pulse [K]:"<< " " << input::T_pulse << "\t"<< std::endl;
    inFile >> s1>> input::pulse_duration;
    std::cout <<"Pulse duration [adim.] //Laser dynamics timescale:" << " " << input::pulse_duration << "\t"<< std::endl;
    inFile >> s1>> input::mx_0;
    std::cout <<"Initial mx coordinate [adim.]:"<< " " << input::mx_0 << "\t"<< std::endl;
    inFile >> s1>> input::my_0;
    std::cout <<"Initial my coordinate [adim.]:"<<" " << input::my_0 << "\t"<< std::endl;
    inFile >> s1>> input::mz_0;
    std::cout <<"Initial mz coordinate [adim.]:" << "" << input::mz_0 << "\t"<< std::endl;
    inFile >> s1>> input::m_vs_T_curve;
    std::cout <<"MT_curve [bool]:"<<"\t" << input::m_vs_T_curve << "\t"<< std::endl;
    inFile >> s1>> input::chipar_vs_T_curve;
    std::cout <<"XparT_curve [bool]:"<<"\t" << input::chipar_vs_T_curve << "\t"<< std::endl;
    inFile >> s1>> input::K_vs_T_curve;
    std::cout <<"KT_curve [bool]:"<<"\t" << input::K_vs_T_curve << "\t"<< std::endl;
    inFile >> s1>> input::equilibrate;
    std::cout <<"Initial equilibration [bool]:"<<"\t" << input::equilibrate << "\t"<< std::endl;
    inFile >> s1>> input::laser_dynamics;
    std::cout <<"Laser dynamics [bool]:"<<"\t" << input::laser_dynamics << "\t"<< std::endl;

    std::cout<<std::endl;
    return 0;

}

int read_material_parameters()
{

 std::string s1;
 std::string input_file="matInput.txt";
 std::ifstream inFile(input_file);


if(!inFile)
    {
        std::cout << std::endl << "Failed to open the file: "<< input_file <<std::endl;
        exit(1);
    }

    std::cout<<"Material parameters: "<<"\n";
   	std::cout<<std::endl;
    inFile >> s1>> input::n_at;
    std::cout <<"Number of atoms per unit cell [adim.]:"<< " " << input::n_at << "\t"<< std::endl;
    inFile >> s1>> input::a;
    std::cout <<"Lattice spacing [m]:"<<  " " << input::a << "\t"<< std::endl;
    inFile >> s1>> input::Tc;
    std::cout <<"Curie Temperature [K]:"<< " " << input::Tc << "\t"<< std::endl;
    inFile >> s1>> input::Ms0_CGS;
    std::cout <<"Saturation Magnetisation CGS [emu/cc]:"<< " " << input::Ms0_CGS << "\t"<< std::endl;
    inFile >> s1>> input::K0_CGS;
    std::cout <<"Magnetocrystalline Anisotropy Constant CGS [erg/cc]:"<< " " << input::K0_CGS << "\t"<< std::endl;

    input::Ms0_SI = input::Ms0_CGS*1e3;
    std::cout<<"Saturation Magnetisation SI [A/m]:"<<" "<<input::Ms0_SI<<"\t"<<std::endl;
    input::K0_SI = input::K0_CGS*1e-1;
    std::cout<<"Magnetocrystalline Anisotropy Constant SI [J/m3]:"<<" "<<input::K0_SI<<"\t"<<std::endl;
    input::volume=input::a*input::a*input::a;
    std::cout<<"Volume [m3]:"<<" "<<input::volume<<"\t"<<std::endl;
    input::mu_s=input::volume*input::Ms0_SI/input::n_at;
    std::cout<<"Atomic magnetic moment [J/T]:"<<" "<<input::mu_s<<"\t"<<std::endl;
	
	inFile >> s1>> input::lambda;
    std::cout <<"Microscopic damping parameter [adim.]:"<< " " << input::lambda << "\t"<< std::endl;    
    inFile >> s1>> input::eps;
    std::cout <<"Spin-wave correction factor [adim.]:"<< " " << input::eps << "\t"<< std::endl;
    inFile >> s1>> input::ex;
    std::cout <<"Easy-axis X coordinate [adim.]:"<< " " << input::ex << "\t"<< std::endl;
    inFile >> s1>> input::ey;
    std::cout <<"Easy-axis Y coordinate [adim.]:"<< " " << input::ey << "\t"<< std::endl;
    inFile >> s1>> input::ez;
    std::cout <<"Easy-axis Z coordinate [adim.]:"<< " " << input::ez << "\t"<< std::endl;
   	std::cout<<std::endl;

    return 0;

}	
    

}


//---End of input.cpp file.

