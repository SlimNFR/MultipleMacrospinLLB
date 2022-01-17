//---This is the output.h file.
#pragma once

#ifndef OUTPUT_H
#define OUTPUT_H

//---Standard libraries
#include<fstream>
#include<iostream>
#include<vector>

//---User-defined libraries


//---Namespace output

namespace output{
	
//---Variables
extern std::vector<std::ofstream> files_Meq_temp_NR;//NR=Newton-Raphson
extern std::vector<std::ofstream> files_Meq_temp_CS;//CS=cubicspline
extern std::vector<std::ofstream> files_X_temp;//X=susceptibility
extern std::vector<std::ofstream> files_K_temp;
extern std::vector<std::ofstream> files_A_temp;
extern std::ofstream file_mx_my_mz_time;
extern std::ofstream file_torques_time;
extern std::ofstream file_energies_cell_time;
extern std::ofstream file_total_energies_time;

//---Functions
int open_files_to_write(int n_materials);
int close_files(int n_materials);
int macrospin_vectors(int n_cells, double time, double T, std::ofstream &f1);
int torques_f(int n_cells, double time, std::ofstream &f1);
int energies_cell_f(int n_cells, double time, std::ofstream &f1);
int total_energies_f(double time, std::ofstream &f1);

	namespace internal{
		int call(int n_cells, double time, double T,
				 std::ofstream &f1,std::ofstream &f2,
				 std::ofstream &f3,std::ofstream &f4);

	}

}

#endif
//---End of output.h file.