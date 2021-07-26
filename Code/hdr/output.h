//---This is the output.h file.
#pragma once

#ifndef OUTPUT_H
#define OUTPUT_H

//---Standard libraries
#include<fstream>
#include<iostream>

//---User-defined libraries


//---Namespace output

namespace output{
	
//---Variables
extern std::ofstream file_Meq_temp_NR;//NR=Newton-Raphson
extern std::ofstream file_Meq_temp_CS;//CS=cubicspline
extern std::ofstream file_X_temp;//X=susceptibility
extern std::ofstream file_K_temp;
extern std::ofstream file_M_time;
//---Functions
int open_files_to_write();
int close_files();

}

#endif
//---End of output.h file.