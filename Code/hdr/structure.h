//---This is the structure.h file.
#pragma once

#ifndef STRUCTURE_H
#define STRUCTURE_H

//---Standard libraries
#include<vector>

//---User-defined libraries
#include"input.h"

//---Namespace macrospin

namespace macrospin{
    
//---Variables
extern std::vector<double>mx, my, mz;


//---Functions
int update_magnetisation(double mx_in,double my_in, double mz_in,
						 double &mx, double &my, double &mz);

}

namespace material{

//---Variables

extern std::vector<double> xcoord;
extern std::vector<double> ycoord;
extern std::vector<double> zcoord;
extern std::vector<int> id;



//---Functions


void generate_f(int n_materials,
                int n_cells,
                std::vector<int>nx, std::vector<int>ny, std::vector<int>nz,
                std::vector<int> &mat_id,
                std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);




}

#endif
//---End of structure.h file.