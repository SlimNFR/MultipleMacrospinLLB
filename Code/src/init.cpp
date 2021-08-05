//---This is the init.cpp file. It defines the init.h functions.


//---Standard libraries
#include<chrono>

//---User-defined libraries
#include"init.h"
#include"tempscaling.h"
#include"input.h"
#include"output.h"
#include"simulation.h"
#include"structure.h"
#include"field.h"
//---Namespace init

namespace init{
	
//---Variables
std::chrono::high_resolution_clock::time_point RUN_TIME_START;
std::chrono::high_resolution_clock::time_point RUN_TIME_END;
double RUN_TIME_TOTAL;


//---Functions
int parameters(int &n_materials)
{	//This function will initialise my simulation parameters

	input::read_material_parameters();
	input::read_simulation_parameters();
	material::generate_crystal_structure_f(input::n_materials,
						 				   input::n_cells, input::nx_cells, input::ny_cells, input::nz_cells,
						 				   material::id, material::xcoord, material::ycoord, material::zcoord);
	material::create_interaction_list(input::n_cells, material::xcoord, material::ycoord, material::zcoord,
									  material::interaction_list, material::start_neighbours, material::end_neighbours);

	
	/*
	std::cout<<"Start neighbours:"<<" ";
	for(int i=0;i<input::n_cells; i++)
	{

		std::cout<<material::start_neighbours[i]<<" ";
	}		
	std::cout<<"\n";
	
	
	
	std::cout<<"End neighbours:"<<" ";
	for(int i=0;i<input::n_cells; i++)
	{

		std::cout<<material::end_neighbours[i]<<" ";
	}		
	std::cout<<"\n";
	
	
	std::cout<<"Interaction list "<<" ";
	for(int i=0;i<material::interaction_list.size(); i++)
	{

		std::cout<<material::interaction_list[i]<<" ";
	}		
	std::cout<<"\n";
	*/

	macrospin::internal::alloc_memory(input::n_cells);
	macrospin::internal::set_initial_config(input::n_cells,
                          		    		input::mx_0, input::my_0, input::mz_0,
                           		  			macrospin::mx_0, macrospin::my_0, macrospin::mz_0,
                           		  			material::id);
	field::internal::alloc_memory(input::n_cells);
	for(int material=0; material<n_materials; material++)
	{
		tempscaling::internal::obtain_interpolation_polynome_mVsT_data(material); //initialise me_vs_T interpolation function	
	}
	
	return 0;
}

int files(int &n_materials)
{// This function will create all the required output files

	output::open_files_to_write(n_materials);

	return 0;
}

int sim(int &n_materials)

{//This function will initialise the chosen simulation

	if(input::m_vs_T_curve==true)
	{
		for(int material=0; material<n_materials; material++)
		{

			tempscaling::internal::call_mVsT_sim(material, output::files_Meq_temp_NR[material], output::files_Meq_temp_CS[material]);

		}
		
	}

	if(input::chipar_vs_T_curve==true)
	{
		for(int material=0; material<n_materials; material++)
		{
			tempscaling::internal::call_chiparVsT_sim(material, output::files_X_temp[material]);
		}
	}

	if(input::K_vs_T_curve==true)
	{
		for(int material=0; material<n_materials; material++)
		{
			tempscaling::internal::call_KVsT_sim(material, output::files_K_temp[material]);
		}
	}

	if(input::A_vs_T_curve==true)
	{
		for(int material=0; material<n_materials; material++)
		{
			tempscaling::internal::call_AvsT_sim(material, output::files_A_temp[material]);
		}
	}

	if(input::equilibrate==true)
	{
		for(int material=0; material<n_materials; material++)
		{
			tempscaling::internal::calc_parameters_at_T(material);
		}
		simulation::equilibrate_system(input::n_cells,
									   input::gamma, input::alpha_par, input::alpha_perp,
									   macrospin::mx_0, macrospin::my_0, macrospin::mz_0,
									   macrospin::mx, macrospin::my, macrospin::mz,
									   field::Bx_eff,field::By_eff, field::Bz_eff,
									   field::torque_mod,
									   material::id,
									   simulation::equil_sim_time, input::T, input::TOL_EQ,
							   		   input::t_min_equil, input::t_max_equil, input::delta_t_equil, input::timescale_equil, 
							   		   output::file_mx_my_mz_time);

	}

	if(input::laser_dynamics == true)
	{

		simulation::squared_pulse_dynamics(input::n_cells,
									       input::gamma, input::alpha_par, input::alpha_perp,
									       macrospin::mx, macrospin::my, macrospin::mz,
									       macrospin::mx, macrospin::my, macrospin::mz,
									       field::Bx_eff,field::By_eff, field::Bz_eff,
									       field::torque_mod,
									   	   material::id,
										   simulation::laser_sim_time, simulation::equil_sim_time, 
										   input::T, input::TOL_LD,
										   input::t_min_laser_dynamics, input::t_max_laser_dynamics,input::delta_t_laser_dynamics, input::timescale_laser_dynamics,
										   input::pulse_duration, input::T_pulse, 
										   output::file_mx_my_mz_time);

	}

	return 0;
	
}



}


namespace init{


	namespace internal{
		int run()
		{//This function will run the code
			init::RUN_TIME_START=std::chrono::high_resolution_clock::now();
			init::parameters(input::n_materials);
			init::files(input::n_materials);
			init::sim(input::n_materials);
			return 0;
		}

		int end()
		{

			output::close_files(input::n_materials);
			init::RUN_TIME_END=std::chrono::high_resolution_clock::now();
			init::RUN_TIME_TOTAL=std::chrono::duration_cast<std::chrono::nanoseconds>(init::RUN_TIME_END - init::RUN_TIME_START).count();
			std::cout<<"Total run time: "<<init::RUN_TIME_TOTAL*1e-9<<" s"<<"\n";
			return 0;

		}	
	}
}




//---End of init.cpp file.

