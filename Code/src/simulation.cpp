//---This is the simulation.cpp file. It defines the simulation.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>

//---User-defined libraries
#include"simulation.h"
#include"solver.h"
#include"equation.h"
#include"field.h"
#include"tempscaling.h"
#include"utils.h"
#include"structure.h"

//---Namespace simulation

namespace simulation{
	
//---Variables

double laser_sim_time = 0.0;
double equil_sim_time = 0.0;
double total_sim_time = laser_sim_time + equil_sim_time;

//---Functions


int force_DW_formation_f(std::vector<double> m_e, std::vector<double> &mx, std::vector<double> &my, std::vector<double> &mz,
						 std::vector<int>material_id)
{	//This function will force the macrospin vectors at the boundaries to be antiparallel
	//I assume the two macrospins are from the same material
	int mat_id;
	int id;
	//double ax, ay, az;
	//double theta = 25.0*M_PI/180.0;
 	//double phi = 0.0*M_PI/180.0;

	
 	

	//ax=mx[id]*cos(theta) + mz[id]*sin(theta);
	//ay=my[id];
 	//az=-mx[id]*sin(theta) + mz[id]*cos(theta);	
	
	id=0;
	mat_id = material_id[id];
	mx[id] = 0.0;
	my[id] = 0.0;
	mz[id] = m_e[mat_id];	

	field::torque_mod[id]=0.0;
	//field::Bx_eff[id]=field::By_eff[id]=field::Bz_eff[id]=0.0;

 	id=mx.size()-1;
	mat_id = material_id[id]; //I assume it's the first and last cell are from the same material
	mx[id] = 0.0;
	my[id] = 0.0;
	mz[id] = -m_e[mat_id];

	field::torque_mod[id]=0.0;
	//field::Bx_eff[id]=field::By_eff[id]=field::Bz_eff[id]=0.0;

	return 0;

}

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
					   	   std::ofstream &f1)
{	
	//This function simulates laser-induced dynamics via a squared pulse

	//Temporary variables
	int t;
	int mat; //this will save the material id
	int max_torque_spin_id;//This will save the id of the spin with the maximum torque acting on it
	double max_torque_mod; //This will hold the maximum torque modulus
	double T_init = T; //save the temperature before laser pulse is applied
	bool ONCE = false;
	
	//Set initial magnetisation coordinates.
	for(int cell = 0; cell<n_cells; cell++)
	{
		mat=material_id[cell];

		mx_n1[cell] = mx_0[cell];
		my_n1[cell] = my_0[cell];
		mz_n1[cell] = mz_0[cell];

	}

		
	//Time loop
	for(t=t_start; t<=t_end; t=t+t_step) //Loop time
	{	
		LD_REAL_t = t*timescale + EQ_REAL_t; //If there was an equilibration time period previously I need to take it into account
			if(t-t_start<=pulse_duration) //For the duration of the pulse..
			{			
				if(t-t_start==0) // Set the medium temperature and calculate the temp.dependent parameters. 
								 //This is done only at the first time instance of the pulse
				{
					T = T_pulse;
					for(int material=0; material<input::n_materials; material++)
					{tempscaling::internal::calc_parameters_at_T(material);}
					
				}

			}
			else
			{
				if(t-t_start>pulse_duration && ONCE == false) // exactly 1 timestep after the pulse ended calculate new temp. dep. params.
				{
					T=T_init;
					for(int material=0; material<input::n_materials; material++)
					{tempscaling::internal::calc_parameters_at_T(material);}
					ONCE = true; //This variable will make sure I calculate the temp. dependent parameters only once after the pulse is taken away. 
								 //The temperature will stay constant so a one-time calculation will be enough.
				}

			}


		field::calculate(); //Compute the field at each new time step
		utils::max_element_1D_vec(torque_mod,max_torque_mod,max_torque_spin_id,false);//Calculate the maximum torque modulus
		

		for(int cell=0; cell<n_cells; cell++) //Loop cells for each time-step
		{	

			mat = material_id[cell]; //Get material id
			f1<<LD_REAL_t<<" "<<mx_n1[cell]<<" "<<my_n1[cell]<<" "<<mz_n1[cell]<<" "<<T<<" "<<"\n";//Print time, magnetisation components and temperature

			//Get the next magnetisation value
			solver::heun_scheme_step(equation::LLB_classic,
									 input::remove_precession_term,
									 mx_0[cell],my_0[cell],mz_0[cell],
									 gamma, alpha_par[mat], alpha_perp[mat],
									 t_step, timescale,
									 Bx_eff[cell],By_eff[cell],Bz_eff[cell],
									 mx_n1[cell], my_n1[cell], mz_n1[cell]);

			//Set the new initial condition
			mx_0[cell] = mx_n1[cell];
			my_0[cell] = my_n1[cell];
			mz_0[cell] = mz_n1[cell];

		}

		if(max_torque_mod<TOL &&  (t-t_start) > pulse_duration )break; //Wait until after the pulse is taken off to apply the breaking condition.
																			  //Otherwise the code might stop just because I start from an equilibrium position.

	}




	return 0;
}

int equilibrate_system(int n_cells,
					   double gamma, std::vector<double> &alpha_par, std::vector<double>&alpha_perp,
					   std::vector<double> mx_0, std::vector<double> my_0, std::vector<double> mz_0,
					   std::vector<double> &mx_n1, std::vector<double> &my_n1, std::vector<double> &mz_n1,
					   std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff,
					   std::vector<double> &torque_mod,
					   std::vector<int> material_id,
					   double &EQ_REAL_t, double T, double TOL,
					   int t_start, int t_end, int t_step, double timescale,
					   std::ofstream &f1)
{	//This function equilibrates the system for a given initial temperature and effective field.
	//The field is calculated as the magnetisation vector changes in time. 

	//Temporary variables
	int mat;//This will save the material id
	int max_torque_spin_id;// This will save the id of the spin with the maximum torque acting on it
	double max_torque_mod;

	//Set initial magnetisation coordinates.
	for(int cell = 0; cell<n_cells; cell++)
	{
		mat=material_id[cell];

		mx_n1[cell] = mx_0[cell];
		my_n1[cell] = my_0[cell];
		mz_n1[cell] = mz_0[cell];

	}
	
	//Time loop
	for(int t=t_start; t<=t_end; t=t+t_step)//Loop time
	{	//std::cout<<"TIME: "<<t<<" TIME_END:"<<t_end<<"\n";

		
		//Force edge spins to be AP		
		if(input::force_DW_formation)simulation::force_DW_formation_f(input::m_e,
																  	  mx_n1, my_n1, mz_n1,
																  	  material_id);	
	
		EQ_REAL_t = t*timescale; //This is the real equilibration time obtained multiplying the imaginary time t by the associated timescale 
		field::calculate(); //Compute the field at each new timestep

		utils::max_element_1D_vec(torque_mod, max_torque_mod,max_torque_spin_id, false); //Calculate the maximum torque in the system

		for(int cell=0; cell<n_cells; cell++)//Loop cells for every timestep
		{
			
			f1<<EQ_REAL_t<<" "<<cell<<" "<<mx_n1[cell]<<" "<<my_n1[cell]<<" "<<mz_n1[cell]<<" "<<T<<" "<<"\n"; //Print time, magnetisation components and temperature

			//Get the next magnetisation value
			solver::heun_scheme_step(equation::LLB_classic,
									 input::remove_precession_term,
									 mx_0[cell],my_0[cell],mz_0[cell],
									 gamma, alpha_par[mat], alpha_perp[mat],
									 t_step, timescale,
									 Bx_eff[cell],By_eff[cell],Bz_eff[cell],
									 mx_n1[cell], my_n1[cell], mz_n1[cell]); 

			//Set new initial condition
			mx_0[cell] = mx_n1[cell];
			my_0[cell] = my_n1[cell];
			mz_0[cell] = mz_n1[cell];

		}
		if(max_torque_mod<TOL)break;//Check for stopping condition

	}
	
	return 0;
}

}


//---End of simulation.cpp file.

