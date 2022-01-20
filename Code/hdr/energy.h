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
extern double anis_total;
extern double exchange_total;
extern double zeeman_total;
extern double total;


//Energies per cell
extern std::vector<double> anis_cell;
extern std::vector<double> exchange_cell;
extern std::vector<double> zeeman_cell;
extern std::vector<double> total_cell;

//---Functions
int anis_cell_f(int n_cells,
                std::vector<double> K, std::vector<double> V,
                std::vector<double> ex, std::vector<double> ey, std::vector<double> ez,
                std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                std::vector<int> mat_id,
                std::vector<double> &anis_cell);


int zeeman_cell_f(int n_cells,
                  double B, double bx, double by, double bz,
                  std::vector<double> Ms, std::vector<double> V,
                  std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                  std::vector<int> mat_id,
                  std::vector<double> &zeeman_cell);


int exchange_cell_f(int n_cells,double lengthscale,
                    std::vector<unsigned long long int> macrocell_size, std::vector<std::vector<double>>A_T_matrix,
                    std::vector<int> int_list, std::vector<int> start_neighbours, std::vector<int> end_neighbours,
                    std::vector<int> material_id,
                    std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
                    std::vector<double>& exc_cell);

 int total_cell_f(int n_cells,
                  std::vector<double> anis_cell,
                  std::vector<double> zeeman_cell,
                  std::vector<double> exchange_cell,
                  std::vector<double> &total_cell);


 int anis_total_f(int n_cells,
                  std::vector<double> anis_cell,
                  double &anis_total);

 int exchange_total_f(int n_cells,
                      std::vector<double> exc_cell,
                      double &exc_total);

 int zeeman_total_f(int n_cells,
                    std::vector<double> zeeman_cell,
                    double &zeeman_total);

 int total_f(double ani,
             double exc,
             double zee,
             double & total);


  namespace internal {

      int alloc_memory(int n_cells);
      int calculate();
  }

}

#endif
//---End of energy.h file.