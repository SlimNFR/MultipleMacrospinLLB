//---This is the energy.h file.
#pragma once

#ifndef ENERGY_H
#define ENERGY_H

//---Standard libraries
#include<vector>

//---User-defined libraries


//---Namespace energy

namespace energy{
	
//---Variables
extern double uniax_anis;
extern double zeeman;

//---Functions
double uniax_anis_f(int n_cells,
					std::vector<double> K, std::vector<double> V,
                    std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                    std::vector<double> ex, std::vector<double> ey, std::vector<double> ez,
                    std::vector<int> mat_id);


double zeeman_f(int n_cells,
                 double B, double bx, double by, double bz,
                 std::vector<double> Ms, std::vector<double> V,
                 std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                 std::vector<double> mat_id);

}

#endif
//---End of energy.h file.