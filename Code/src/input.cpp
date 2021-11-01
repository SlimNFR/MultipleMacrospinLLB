//---This is the input.cpp file. It defines the input.h functions.


//---Standard libraries
#include<cmath>
#include<iostream>
#include<fstream>
#include<string>
#include<deque>

//---User-defined libraries
#include"input.h"

//---Namespace input

namespace input{
	

//---Variables

// //---Global constants


// //---Material parameters
int n_materials; //No. of materials in my system
double lengthscale; //The lengthscale can be nm, microns, mm and so on..
std::vector<unsigned long long int>length; //Material lengths
std::vector<unsigned long long int>width; //Material widths
std::vector<unsigned long long int>height; //Material heights
std::vector<unsigned long long int>unitcell_size; //This sets the a dimension of the unitcell
std::vector<double>unitcell_volume; //This sets the volume of the unitcell
std::vector<unsigned long long int>macrocell_size; //Assuming a cubic cell type, this sets the "a" dimension of the macrocell
std::vector<double>macrocell_volume; // Macrocell volume: [m^3]
std::vector<unsigned long long int>material_total_cells; //number of total cells per material
std::vector<unsigned long long int>material_nx_cells; //Total cells on x per material
std::vector<unsigned long long int>material_ny_cells; //Total cells on y per material
std::vector<unsigned long long int>material_nz_cells; //Total cells on z per material
int nx_cells;//Total cells on x
int ny_cells;//Total cells on y
int nz_cells;//Total cells on z
int n_cells; //Total number of cells in the structure
int n_at; //Number of atoms per unit cell: adim.

std::vector<double>Tc; // Curie temperature: [K]
std::vector<double>Ms0_CGS; // Saturation magnetisation: [emu/cc]
std::vector<double>K0_CGS; // Magnetocrystlline first-order anis. constant: [erg/cc]
std::vector<double>Ms0_SI; //Saturation magnetisation: [A/m]
std::vector<double>K0_SI; //: Magnetocrystlline first-order anis. constant: [J/m^3]
std::vector<double> A0_list; //0K Exchange stiffness list [J/m]
std::vector<std::vector<double>> A0_matrix;//0K Exchange stiffness matrix [J/m]
std::vector<std::vector<double>> A_T_matrix;//Temperature dependent exchange stiffness matrix [J/m]
std::vector<double>m_e; //Equilibrium magnetisation: []
std::vector<double>Ms_T; //Saturation magnetisation at T
std::vector<double>K_T; //Magnetocrystalline anisotropy constant at T
std::vector<double>mu_s; // Atomic magnetic moment: [J/T]
std::vector<double>lambda; // Microscopic coupling constant: adimensional
std::vector<double>chi_par; //Parallel susceptibility: []
std::vector<double>chi_perp; //Perpendicular susceptibility: []
std::vector<double>eps; //SW correction factor: adim.
std::vector<double>ex; // H_ani_x
std::vector<double>ey; // H_ani_y
std::vector<double>ez; // H_ani_z
//---Initial conditions
std::vector<double>mx_0;
std::vector<double>my_0;
std::vector<double>mz_0;
bool read_config_from_file;

// //---Simulation paramters
double T; // Magnetic moment temperature (electron temperature in this model): [K]
double B_app; // External field amplitude: [T]
double B_theta;//angle with respect to EA
double B_phi; //Angle between in-plane field projection and the Ox axis
double bx; // H_ext_x/H_ext 
double by; // H_ext_y/H_ext
double bz; // H_ext_z/H_ext
std::vector<double>alpha_par; //Longitudinal damping parameter: adimensional. Avail. only for T<Tc
std::vector<double>alpha_perp; //Transversal damping parameter: adimensional. Avail. only for T<Tc

int t_min_equil; // Time will be given as an integer. It needs to be multiplied with timescale_equil to obtain the real-time.
int t_max_equil;
int t_step_output_equil;
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



//---Simulation

bool m_vs_T_curve;
bool chipar_vs_T_curve;
bool K_vs_T_curve;
bool A_vs_T_curve;
bool equilibrate;
bool laser_dynamics;
bool force_DW_formation;
bool remove_precession_term;

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
    std::cout <<"Applied field theta angle [degr.]:"<< " " << input::B_theta << "\t"<< std::endl;
    input::B_theta=input::B_theta*M_PI/180.0;

    inFile >> s1>> input::B_phi;
    std::cout <<"Applied field phi angle [degr.]:"<< " " << input::B_phi << "\t"<< std::endl;  
    input::B_phi=input::B_phi*M_PI/180.0;

    input::bx= sin(input::B_theta)*cos(input::B_phi); // H_ext_x/H_ext 
	input::by = sin(input::B_theta)*sin(input::B_phi); // H_ext_y/H_ext
	input::bz = cos(input::B_theta); // H_ext_z/H_ext
	std::cout<<"Applied field X coordinate [adim.]:"<<" "<<input::bx<<"\n";
	std::cout<<"Applied field Y coordinate [adim.]:"<<" "<<input::by<<"\n";
	std::cout<<"Applied field Z coordinate [adim.]:"<<" "<<input::bz<<"\n";

    input::alpha_par.resize(n_materials);
    input::alpha_perp.resize(n_materials);

  	inFile >> s1>> input::t_min_equil;
    std::cout <<"Initial equilibration time start [adim.]:"<< " " << input::t_min_equil << "\t "<< std::endl;
    inFile >> s1>> input::t_max_equil;
    std::cout <<"Initial equilibration time stop [adim.]:"<< " " << input::t_max_equil << "\t "<< std::endl;
    inFile >> s1>> input::t_step_output_equil;
    std::cout <<"Time step equilibration output [adim.]:"<< " " << input::t_step_output_equil << "\t "<< std::endl;
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
    
    inFile >> s1>> input::m_vs_T_curve;
    std::cout <<"MT_curve [bool]:"<<"\t" << input::m_vs_T_curve << "\t"<< std::endl;
    inFile >> s1>> input::chipar_vs_T_curve;
    std::cout <<"XparT_curve [bool]:"<<"\t" << input::chipar_vs_T_curve << "\t"<< std::endl;
    inFile >> s1>> input::K_vs_T_curve;
    std::cout <<"KT_curve [bool]:"<<"\t" << input::K_vs_T_curve << "\t"<< std::endl;
    inFile >> s1>> input::A_vs_T_curve;
    std::cout <<"AT_curve [bool]:"<<"\t" << input::A_vs_T_curve << "\t"<< std::endl;
    inFile >> s1>> input::equilibrate;
    std::cout <<"Initial equilibration [bool]:"<<"\t" << input::equilibrate << "\t"<< std::endl;
    inFile >> s1>> input::laser_dynamics;
    std::cout <<"Laser dynamics [bool]:"<<"\t" << input::laser_dynamics << "\t"<< std::endl;
    inFile >> s1>> input::force_DW_formation;
    std::cout <<"Force DW formation [bool]:"<<"\t" <<input::force_DW_formation<< "\t"<< std::endl;
    inFile >> s1>> input::remove_precession_term;
    std::cout <<"Remove precession term from LLB [bool]:"<<"\t" <<input::remove_precession_term<< "\t"<< std::endl;


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

    inFile >> s1>> input::n_materials;
    std::cout <<"Number of materials simulated [adim.]:"<< " " << input::n_materials << "\t"<< std::endl;

    inFile >> s1>> input::lengthscale;
    std::cout <<"Lengthscale [m]:"<< " " << input::lengthscale << "\t"<< std::endl;
    
    inFile>>s1;
    input::length.resize(n_materials);
    std::cout<<"Material lengths [m]:"<<" ";
    for (int i = 0; i <input::n_materials;i++)
    {
        inFile>>input::length[i];
        std::cout<<input::length[i]*input::lengthscale<<" ";
    }
    std::cout<<"\n";

    inFile>>s1;
    input::width.resize(n_materials);
    std::cout<<"Material widths [m]:"<<" ";
    for (int i = 0; i <input::n_materials;i++)
    {
        inFile>>input::width[i];
        std::cout<<input::width[i]*input::lengthscale<<" ";
    }
    std::cout<<"\n";

    inFile>>s1;
    input::height.resize(n_materials);
    std::cout<<"Material heights [m]:"<<" ";
    for (int i = 0; i <input::n_materials;i++)
    {
        inFile>>input::height[i];
        std::cout<<input::height[i]*input::lengthscale<<" ";
    }
    std::cout<<"\n";

    inFile>>s1;
    input::unitcell_size.resize(n_materials,0);
    input::unitcell_volume.resize(n_materials,0);
    input::macrocell_size.resize(n_materials,0);
    input::macrocell_volume.resize(n_materials,0);
    input::material_total_cells.resize(n_materials,0);
    input::material_nx_cells.resize(n_materials,0);
    input::material_ny_cells.resize(n_materials,0);
    input::material_nz_cells.resize(n_materials,0);
    input::nx_cells=0;
    input::ny_cells=0;
    input::nz_cells=0;
    input::n_cells=0;



    std::cout<<"Unitcell size/a dimension [m]:"<<" ";
    for (int i = 0; i <input::n_materials;i++)
    {
        inFile>>input::unitcell_size[i];
        std::cout<<(double)input::unitcell_size[i]*input::lengthscale<<" ";


        //At the same time, I can calculate the unitcellvolume
        input::unitcell_volume[i]=pow(input::unitcell_size[i]*input::lengthscale,3.0);
        
    }
    std::cout<<"\n";

    std::cout<<"Unitcell volume [m3]:"<<" ";
    for (int i = 0; i < input::n_materials; i++)
    {
        std::cout<<input::unitcell_volume[i]<<" ";
    }
    std::cout<<"\n";

    inFile>>s1;
    std::cout<<"Macrocell size/a dimension [m]:"<<" ";
    for (int i = 0; i <input::n_materials;i++)
    {
        inFile>>input::macrocell_size[i];
        std::cout<<(double)input::macrocell_size[i]*input::lengthscale<<" ";


        //At the same time, I can calculate the macrocell volume for each material and the number of total cells in each material
       input::macrocell_volume[i]=pow(input::macrocell_size[i]*input::lengthscale,3.0);
       input::material_total_cells[i]=(input::length[i]*input::width[i]*input::height[i])/(input::macrocell_size[i]*input::macrocell_size[i]*input::macrocell_size[i]);
       std::cout<<"MATERIAL_TOTAL_CELLS: "<<input::material_total_cells[i]<<" |LENGTH: "<<input::length[i]<<" |WIDTH: "<<input::width[i]
                <<" |HEIGHT: "<<input::height[i]<<" |MACROCELL_SIZE:"<<input::macrocell_size[i]<<"\n";
        std::cout<<"TOTAL VOLUME: "<<input::length[i]*input::width[i]*input::height[i]<<"\n";
        std::cout<<"macrocell VOLUME: "<<input::macrocell_size[i]*input::macrocell_size[i]*input::macrocell_size[i]<<"\n";

        //I also calculate the total cells on x,y,z per material

        input::material_nx_cells[i]= input::length[i]/input::macrocell_size[i];
        input::material_ny_cells[i]= input::width[i]/input::macrocell_size[i];
        input::material_nz_cells[i]= input::height[i]/input::macrocell_size[i];

        //The total number of cells in the system
        input::nx_cells +=input::material_nx_cells[i];
        input::ny_cells +=input::material_ny_cells[i];
        input::nz_cells +=input::material_nz_cells[i];
        input::n_cells += input::material_total_cells[i]; 
        std::cout<<"N_CELLS: "<<n_cells<<"\n";

    }
    
    std::cout<<"\n";


    std::cout<<"Macrocell volume [m3]:"<<" ";
    for (int i = 0; i < input::n_materials; i++)
    {
        std::cout<<input::macrocell_volume[i]<<" ";
    }
    std::cout<<"\n";

    
    std::cout<<"Total number of cells [adim.]:"<<input::n_cells<<"\n";

    std::cout<<"No. of cells per material [adim.]:"<<" ";
    for (int i = 0; i < input::n_materials; i++)
    {
        std::cout<<input::material_total_cells[i]<<" ";
    }
    std::cout<<"\n";

    std::cout<<"No. of cells per material on x [adim.]:"<<" ";
    for (int i = 0; i < input::n_materials; i++)
    {
        std::cout<<input::material_nx_cells[i]<<" ";
    }
    std::cout<<"\n";

    std::cout<<"No. of cells per material on y [adim.]:"<<" ";
    for (int i = 0; i < input::n_materials; i++)
    {
        std::cout<<input::material_ny_cells[i]<<" ";
    }
    std::cout<<"\n";

    std::cout<<"No. of cells per material on z [adim.]:"<<" ";
    for (int i = 0; i < input::n_materials; i++)
    {
        std::cout<<input::material_nz_cells[i]<<" ";
    }
    std::cout<<"\n";

    


    inFile >> s1>> input::n_at;
    std::cout <<"Number of atoms per unit cell [adim.]:"<< " " << input::n_at << "\t"<< std::endl;

    
    inFile >> s1;
    input::Tc.resize(n_materials);
    std::cout <<"Curie Temperature [K]:"<< " " ;
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>>input::Tc[i];
        std::cout<< input::Tc[i] << " ";
    }
    std::cout<<"\n";
    
    inFile >> s1;
    input::Ms0_CGS.resize(n_materials);
    std::cout <<"Saturation Magnetisation CGS [emu/cc]:"<< " ";
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>> input::Ms0_CGS[i];
        std::cout<< input::Ms0_CGS[i] << " ";
    }
    std::cout<<"\n";


    inFile >> s1;
    input::K0_CGS.resize(n_materials);
    std::cout <<"Magnetocrystalline Anisotropy Constant CGS [erg/cc]:"<<" ";
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>> input::K0_CGS[i];
        std::cout<< input::K0_CGS[i] << " ";
    }
    std::cout<<"\n";
    

    input::Ms0_SI.resize(n_materials);
    std::cout<<"Saturation Magnetisation SI [A/m]:"<<" ";
    for(int i=0; i<input::n_materials; i++)
    {
        input::Ms0_SI[i] = input::Ms0_CGS[i]*1e3;
        std::cout<< input::Ms0_SI[i] << " ";
    }
    std::cout<<"\n";

    

    input::K0_SI.resize(n_materials);
    std::cout<<"Magnetocrystalline Anisotropy Constant SI [J/m3]:"<<" ";
    for(int i=0; i<input::n_materials; i++)
    {
        input::K0_SI[i] = input::K0_CGS[i]*1e-1;
        std::cout<< input::K0_SI[i] << " ";
    }
    std::cout<<"\n";

    inFile>>s1;
    input::A0_list.resize(2*n_materials-1);
    std::cout<<"Exchange stiffness list [J/m]:"<<" ";
    for(int i=0; i<2*n_materials-1; i++)
    {
        inFile>>input::A0_list[i];
        std::cout<<input::A0_list[i]<<" ";

    }
    std::cout<<"\n";

    input::generate_exchange_matrix(input::n_materials,
                                    input::A0_list,
                                    input::A0_matrix);

    input::A_T_matrix.resize(n_materials, std::vector<double>(n_materials));

    std::cout<<"Exchange stiffness matrix [J/m]:"<<"\n";
    for(int i=0; i<n_materials; i++)
    {
        for (int j=0; j<n_materials;j++)
        {
            std::cout<<A0_matrix[i][j]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";

    input::m_e.resize(n_materials);
    input::Ms_T.resize(n_materials);
    input::K_T.resize(n_materials);

    input::mu_s.resize(n_materials);
    std::cout<<"Atomic magnetic moment [J/T]:"<<" ";
    for(int i=0; i<input::n_materials; i++)
    {
        input::mu_s[i] = input::unitcell_volume[i]*input::Ms0_SI[i]/(double)input::n_at;
        std::cout<< input::mu_s[i] << " ";
    }
    std::cout<<"\n";    

    inFile >> s1;
    input::lambda.resize(n_materials);
    std::cout <<"Microscopic damping parameter [adim.]:"<< " ";
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>> input::lambda[i];
        std::cout<< input::lambda[i] << " ";
    }
    std::cout<<"\n";  


    input::chi_par.resize(n_materials);
    input::chi_perp.resize(n_materials);

    inFile >> s1;
    input::eps.resize(n_materials);
    std::cout <<"Spin-wave correction factor [adim.]:"<< " ";
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>>input::eps[i];
        std::cout<< input::eps[i] << " ";
    }
    std::cout<<"\n";  	

    inFile >> s1;
    input::ex.resize(n_materials);
    std::cout <<"Easy-axis X coordinate [adim.]:"<< " ";
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>>input::ex[i];
        std::cout<< input::ex[i] << " ";
    }
    std::cout<<"\n";    

    inFile >> s1;
    input::ey.resize(n_materials);
    std::cout <<"Easy-axis Y coordinate [adim.]:"<< " ";
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>>input::ey[i];
        std::cout<< input::ey[i] << " ";
    }
    std::cout<<"\n";


    inFile >> s1;
    input::ez.resize(n_materials);
    std::cout <<"Easy-axis Z coordinate [adim.]:"<< " ";
    for(int i=0; i<input::n_materials; i++)
    {
        inFile>>input::ez[i];
        std::cout<< input::ez[i] << " ";
    }
    std::cout<<"\n";

    inFile >> s1;
    std::cout <<"Initial mx coordinate [adim.]:"<< " ";
    input::mx_0.resize(input::n_materials);
    for(int i=0; i<n_materials;i++)
    {
        inFile>>input::mx_0[i];
        std::cout<< input::mx_0[i] << " ";
    }
    std::cout<<"\n";

    inFile >> s1;
    std::cout <<"Initial my coordinate [adim.]:"<< " ";
    input::my_0.resize(input::n_materials);
    for(int i=0; i<n_materials;i++)
    {
        inFile>>input::my_0[i];
        std::cout<< input::my_0[i] << " ";
    }
    std::cout<<"\n";

    inFile >> s1;
    std::cout <<"Initial mz coordinate [adim.]:"<< " ";
    input::mz_0.resize(input::n_materials);
    for(int i=0; i<n_materials;i++)
    {
        inFile>>input::mz_0[i];
        std::cout<< input::mz_0[i] << " ";
    }

    std::cout<<"\n";
    inFile >> s1>> input::read_config_from_file;
    std::cout <<"Read spin configuration from file [bool]:"<<"\t" <<input::read_config_from_file<< "\t"<< std::endl;
    
    std::cout<<"\n";

    return 0;

}	



 int generate_exchange_matrix(int n_materials,
                              std::vector<double>A0_list,
                              std::vector<std::vector<double> >&A0_matrix)
 {// This function generates the exchange stiffness matrix from the exchange stiffness list given as input
  A0_matrix.resize(n_materials, std::vector<double>(n_materials));

  for(int i=0;i<n_materials;i++) //Loop each material
  {
   for(int j=0;j<n_materials;j++) //For every other material
   {
    if(abs(j-i)>1.5)//This determines whether the two materials interact
    {//Currently, only two succesive materials interact. 
     A0_matrix[i][j]=0.0;
     continue;
    }
    else
    {
     A0_matrix[i][j] = A0_list[i+j];
    }
   }
  }

  return 0;
 }
    

}


//---End of input.cpp file.

