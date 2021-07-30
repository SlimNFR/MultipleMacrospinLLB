//---This is the output.cpp file. It defines the output.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>
#include<string>

//---User-defined libraries
#include"output.h"
#include"input.h"

//---Namespace output

namespace output{
	
//---Variables
std::vector<std::ofstream> files_Meq_temp_NR;//NR=Newton-Raphson
std::vector<std::ofstream> files_Meq_temp_CS;//CS=cubicspline
std::vector<std::ofstream> files_X_temp;//X=susceptibility
std::vector<std::ofstream> files_K_temp;
std::ofstream file_mx_my_mz_time;

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

	if(input::equilibrate == true || input::laser_dynamics == true )
	{

		output::file_mx_my_mz_time.open("output_mx_my_mz_time.txt", std::ofstream::out);	

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


	if(input::equilibrate == true || input::laser_dynamics == true)
	{

		output::file_mx_my_mz_time.close();

	}


	return 0;

}

}


//---End of output.cpp file.

