//---This is the simulation.cpp file. It defines the simulation.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>
#include<cmath>

//---User-defined libraries
#include"simulation.h"
#include"solver.h"
#include"equation.h"
#include"field.h"
#include"input.h"
#include"particle.h"
#include"tempscaling.h"

//---Namespace simulation

namespace simulation{
	
//---Variables

double laser_sim_time = 0.0;
double equil_sim_time = 0.0;
double total_sim_time = laser_sim_time + equil_sim_time;

//---Functions

int squared_pulse_dynamics(double mx_0, double my_0, double mz_0,
					   	   double &mx_n1, double &my_n1, double &mz_n1,
					   	   double &LD_REAL_t, double EQ_REAL_t,
					   	   double &T, double TOL,
					   	   int t_start, int t_end, int t_step, double timescale,
					   	   int pulse_duration, double T_pulse,
					   	   std::ofstream &f1)
{	
	//This function simulates laser-induced dynamics via a squared pulse

	//Temporary variables
	int t;
	double Bx_eff, By_eff, Bz_eff;
	double torque_mod;
	double gamma, alpha_par, alpha_perp;
	double T_init = T; //save the temperature before laser pulse is applied
	bool ONCE = false;
	
	//Set initial magnetisation coordinates.
	mx_n1 = mx_0;
	my_n1 = my_0;
	mz_n1 = mz_0;

	//Define LLB parameters
	gamma=input::gamma; alpha_par=input::alpha_par; alpha_perp=input::alpha_perp;

	//Time loop
	for(t=t_start; t<=t_end; t=t+t_step) 
	{	
		LD_REAL_t = t*timescale + EQ_REAL_t; //If there was an equilibration time period previously I need to take it into account
		if(t-t_start<=pulse_duration) //For the duration of the pulse..
		{			
			if(t-t_start==0) // Set the medium temperature and calculate the temp.dependent parameters. 
							 //This is done only at the first time instance of the pulse
			{
				T = T_pulse;
				tempscaling::internal::calc_parameters_at_T();
			}

		}
		else
		{
			if(t-t_start>pulse_duration && ONCE == false) // exactly 1 timestep after the pulse ended calculate new temp. dep. params.
			{
				T=T_init;
				tempscaling::internal::calc_parameters_at_T();
				ONCE = true; //This variable will make sure I calculate the temp. dependent parameters only once after the pulse is taken away. 
							 //The temperature will stay constant so a one-time calculation will be enough.
			}

		}

		f1<<LD_REAL_t<<" "<<mx_n1<<" "<<my_n1<<" "<<mz_n1<<" "<<T<<" "<<"\n";//Print time, magnetisation components and temperature

		field::calculate(); //Compute the field at each new step
		//Assign the newly obtained field and torque values to the temporary variables below.
		Bx_eff = field::Bx_eff; By_eff = field::By_eff; Bz_eff = field::Bz_eff;
		torque_mod=field::torque_mod;

		if( fabs(torque_mod)<TOL &&  (t-t_start) > pulse_duration )break; //Wait until after the pulse is taken off to apply the breaking condition.
																		  //Otherwise the code might stop just because I start from an equilibrium position.

		//Get the next magnetisation value
		solver::heun_scheme_step(equation::LLB_classic,
								 mx_0,my_0,mz_0,
								 gamma, alpha_par, alpha_perp,
								 t_step, timescale,
								 Bx_eff,By_eff,Bz_eff,
								 mx_n1, my_n1, mz_n1);

		//Set the new initial condition
		mx_0 = mx_n1;
		my_0 = my_n1;
		mz_0 = mz_n1;


	}




	return 0;
}

int equilibrate_system(double mx_0, double my_0, double mz_0,
					   double &mx_n1, double &my_n1, double &mz_n1,
					   double &EQ_REAL_t, double T, double TOL,
					   int t_start, int t_end, int t_step, double timescale,
					   std::ofstream &f1)
{	//This function equilibrates the system for a given initial temperature and effective field.
	//The field is calculated as the magnetisation vector changes in time. 

	//Temporary variables
	double Bx_eff, By_eff, Bz_eff;
	double torque_mod;
	double gamma, alpha_par, alpha_perp;
	
	//Set initial magnetisation coordinates.
	mx_n1 = mx_0;
	my_n1 = my_0;
	mz_n1 = mz_0;

	//Define LLB parameters
	gamma=input::gamma; alpha_par=input::alpha_par; alpha_perp=input::alpha_perp;

	//Time loop
	for(int t=t_start; t<=t_end; t=t+t_step)
	{	
		EQ_REAL_t = t*timescale; //This is the real equilibration time obtained multiplying the imaginary time t by the associated timescale 
		f1<<EQ_REAL_t<<" "<<mx_n1<<" "<<my_n1<<" "<<mz_n1<<" "<<T<<" "<<"\n"; //Print time, magnetisation components and temperature
		
		
		field::calculate(); //Compute the field at each new step
		//Assign the newly obtained field and torque values to the temporary variables below.
		Bx_eff = field::Bx_eff; By_eff = field::By_eff; Bz_eff = field::Bz_eff;
		torque_mod=field::torque_mod; 

		if(fabs(torque_mod)<TOL)break;//Stopping condition.

		//Get the next magnetisation value
		solver::heun_scheme_step(equation::LLB_classic,
								 mx_0,my_0,mz_0,
								 gamma, alpha_par, alpha_perp,
								 t_step, timescale,
								 Bx_eff,By_eff,Bz_eff,
								 mx_n1, my_n1, mz_n1); 

		//Set new initial condition
		mx_0 = mx_n1;
		my_0 = my_n1;
		mz_0 = mz_n1;


	}
	
	return 0;
}

}


//---End of simulation.cpp file.

