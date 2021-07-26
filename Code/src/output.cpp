//---This is the output.cpp file. It defines the output.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>

//---User-defined libraries
#include"output.h"
#include"input.h"

//---Namespace output

namespace output{
	
//---Variables
std::ofstream file_Meq_temp_NR;//NR=Newton-Raphson
std::ofstream file_Meq_temp_CS;//CS=cubicspline
std::ofstream file_X_temp;//X=susceptibility
std::ofstream file_K_temp;
std::ofstream file_M_time;   

//---Functions
int open_files_to_write()
{

	if(input::m_vs_T_curve == true)
	{
		output::file_Meq_temp_NR.open("output_Meq_temp_NR.txt", std::ofstream::out);
		output::file_Meq_temp_CS.open("output_Meq_temp_CS.txt", std::ofstream::out);
	}

	if(input::chipar_vs_T_curve == true)
	{
		output::file_X_temp.open("output_X_temp.txt", std::ofstream::out);
	}

	if(input::K_vs_T_curve == true)
	{
		output::file_K_temp.open("output_K_temp.txt", std::ofstream::out);
	}

	if(input::equilibrate == true || input::laser_dynamics == true )
	{

		output::file_M_time.open("output_M_time.txt", std::ofstream::out);	

	}

	return 0;
}

int close_files()
{

	if(input::m_vs_T_curve == true)
	{
		output::file_Meq_temp_NR.close();
		output::file_Meq_temp_CS.close();
	}

	if(input::chipar_vs_T_curve == true)
	{
		output::file_X_temp.close();
	}

	if(input::K_vs_T_curve == true)
	{
		output::file_K_temp.close();
	}


	if(input::equilibrate == true || input::laser_dynamics == true)
	{

		output::file_M_time.close();

	}


	return 0;

}

}


//---End of output.cpp file.

