//---This is the init.cpp file. It defines the init.h functions.


//---Standard libraries
#include<chrono>

//---User-defined libraries
#include"init.h"
#include"tempscaling.h"
#include"input.h"
#include"output.h"
#include"simulation.h"
#include"particle.h"

//---Namespace init

namespace init{
	
//---Variables
std::chrono::high_resolution_clock::time_point RUN_TIME_START;
std::chrono::high_resolution_clock::time_point RUN_TIME_END;
double RUN_TIME_TOTAL;


//---Functions
int parameters()
{	//This function will initialise my simulation parameters

	input::read_material_parameters();
	input::read_simulation_parameters();
	tempscaling::internal::obtain_interpolation_polynome_mVsT_data(); //initialise me_vs_T interpolation function
	return 0;
}

int files()
{// This function will create all the required output files

	output::open_files_to_write();

	return 0;
}

int sim()
{//This function will initialise the chosen simulation

	if(input::m_vs_T_curve==true)
	{
		tempscaling::internal::call_mVsT_sim();
	}

	if(input::chipar_vs_T_curve==true)
	{
		tempscaling::internal::call_chiparVsT_sim();
	}

	if(input::K_vs_T_curve==true)
	{

		tempscaling::internal::call_KVsT_sim();
	}

	if(input::equilibrate==true)
	{
		tempscaling::internal::calc_parameters_at_T();
		simulation::equilibrate_system(input::mx_0, input::my_0, input::mz_0,
									   particle::mx, particle::my, particle::mz,
									   simulation::equil_sim_time, input::T, input::TOL_EQ,
							   		   input::t_min_equil, input::t_max_equil, input::delta_t_equil, input::timescale_equil, 
							   		   output::file_M_time);

	}

	if(input::laser_dynamics == true)
	{

		simulation::squared_pulse_dynamics(particle::mx, particle::my, particle::mz,
										   particle::mx, particle::my, particle::mz,
										   simulation::laser_sim_time, simulation::equil_sim_time, 
										   input::T, input::TOL_LD,
										   input::t_min_laser_dynamics, input::t_max_laser_dynamics,input::delta_t_laser_dynamics, input::timescale_laser_dynamics,
										   input::pulse_duration, input::T_pulse, 
										   output::file_M_time);

	}

	return 0;
	
}



}


namespace init{


	namespace internal{
		int run()
		{//This function will run the code
			init::RUN_TIME_START=std::chrono::high_resolution_clock::now();
			init::parameters();
			init::files();
			init::sim();
			return 0;
		}

		int end()
		{

			output::close_files();
			init::RUN_TIME_END=std::chrono::high_resolution_clock::now();
			init::RUN_TIME_TOTAL=std::chrono::duration_cast<std::chrono::nanoseconds>(init::RUN_TIME_END - init::RUN_TIME_START).count();
			std::cout<<"Total run time: "<<init::RUN_TIME_TOTAL*1e-9<<" s"<<"\n";
			return 0;

		}	
	}
}




//---End of init.cpp file.

