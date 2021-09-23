//---This is the output.cpp file. It defines the output.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>
#include<string>

//---User-defined libraries
#include"field.h"
#include"output.h"
#include"input.h"

//---Namespace output

namespace output{
	
//---Variables
std::vector<std::ofstream> files_Meq_temp_NR;//NR=Newton-Raphson
std::vector<std::ofstream> files_Meq_temp_CS;//CS=cubicspline
std::vector<std::ofstream> files_X_temp;//X=susceptibility
std::vector<std::ofstream> files_K_temp;
std::vector<std::ofstream> files_A_temp;
std::ofstream file_mx_my_mz_time;
std::ofstream file_torques_time;

//---Functions
int open_files_to_write(int n_materials)
{

	std::string s1;
	if(input::m_vs_T_curve == true)
	{
		for(int material = 0; material< n_materials; material++)
		{	//see https://stackoverflow.com/questions/31885704/why-can-i-have-a-stdvectorstdofstream-but-not-a-stdvectorstdofstream for explanation
			// on usageof pushback(std::move(temp..));
			s1=std::to_string(material);
			std::ofstream tempNR("output_Meq_temp_NR_mat"+s1+".txt");
			output::files_Meq_temp_NR.push_back(std::move(tempNR));
			std::ofstream tempCS("output_Meq_temp_CS_mat"+s1+".txt");
			output::files_Meq_temp_CS.push_back(std::move(tempCS));

//			output::files_Meq_temp_NR[material].open("output_Meq_temp_NR"+s1+".txt", std::ofstream::out);
//			output::files_Meq_temp_CS[material].open("output_Meq_temp_CS"+s1+".txt", std::ofstream::out);

		}
		
	}

	if(input::chipar_vs_T_curve == true)
	{
		for(int material = 0; material<n_materials; material++)
		{
			s1=std::to_string(material);
			std::ofstream tempX("output_X_temp_mat"+s1+".txt");
			output::files_X_temp.push_back(std::move(tempX));
			//output::files_X_temp[material].open("output_X_temp"+s1+".txt", std::ofstream::out);
		}
		
	}

	if(input::K_vs_T_curve == true)
	{
		for(int material=0; material<n_materials; material++)
		{
			s1=std::to_string(material);
			std::ofstream tempK("output_K_temp_mat"+s1+".txt");
			output::files_K_temp.push_back(std::move(tempK));
			//output::files_K_temp[material].open("output_K_temp"+s1+".txt", std::ofstream::out);	
		}
		
	}

	if(input::A_vs_T_curve == true)
	{
		for(int material=0; material<n_materials; material++)
		{
			s1=std::to_string(material);
			std::ofstream tempA("output_A_temp_mat"+s1+".txt");
			output::files_A_temp.push_back(std::move(tempA));
			//output::files_K_temp[material].open("output_K_temp"+s1+".txt", std::ofstream::out);	
		}
		
	}

	if(input::equilibrate == true || input::laser_dynamics == true )
	{

		output::file_mx_my_mz_time.open("output_mx_my_mz_time.txt", std::ofstream::out);
		output::file_torques_time.open("output_torques_time.txt", std::ofstream::out);	

	}

	return 0;
}

int close_files(int n_materials)
{

	if(input::m_vs_T_curve == true)
	{
		for(int material = 0; material< n_materials; material++)
		{
			output::files_Meq_temp_NR[material].close();
			output::files_Meq_temp_CS[material].close();

		}
	}

	if(input::chipar_vs_T_curve == true)
	{
		for(int material = 0; material<n_materials; material++)
		{
			output::files_X_temp[material].close();
		}
	}

	if(input::K_vs_T_curve == true)
	{
		for(int material=0; material<n_materials; material++)
		{
			output::files_K_temp[material].close();	
		}
	}

	if(input::A_vs_T_curve == true)
	{
		for(int material=0; material<n_materials; material++)
		{
			output::files_A_temp[material].close();	
		}
	}


	if(input::equilibrate == true || input::laser_dynamics == true)
	{

		output::file_mx_my_mz_time.close();
		output::file_torques_time.close();

	}


	return 0;

}

int torques_f(int cell, double time, std::ofstream &f1)
{


	f1<<time<<" "<<cell<<" "<<field::torque_app_x[cell]<<" "<<field::torque_app_y[cell]<<" "<<field::torque_app_z[cell]<<" "<<field::torque_app_mod[cell]<<" "
										 <<field::torque_ani_x[cell]<<" "<<field::torque_ani_y[cell]<<" "<<field::torque_ani_z[cell]<<" "<<field::torque_ani_mod[cell]<<" "
										 <<field::torque_exc_x[cell]<<" "<<field::torque_exc_y[cell]<<" "<<field::torque_exc_z[cell]<<" "<<field::torque_exc_mod[cell]<<" "
										 <<field::torque_lon_x[cell]<<" "<<field::torque_lon_y[cell]<<" "<<field::torque_lon_z[cell]<<" "<<field::torque_lon_mod[cell]<<" "
										 <<field::torque_x[cell]<<" "<<field::torque_y[cell]<<" "<<field::torque_z[cell]<<" "<<field::torque_mod[cell]<<"\n";


	return 0;	
}

}


//---End of output.cpp file.

