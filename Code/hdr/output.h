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
extern std::ofstream file_mx_my_mz_time;
//---Functions
int open_files_to_write(int n_materials);
int close_files(int n_materials);

}

#endif
//---End of output.h file.