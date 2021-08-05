//---This is the field.h file.
#pragma once

#ifndef FIELD_H
#define FIELD_H

//---Standard libraries
#include<vector>

//---User-defined libraries


//---Namespace field

namespace field{
	
//---Variables
extern std::vector<double> Bx_app, By_app, Bz_app;
extern std::vector<double> Bx_ani, By_ani, Bz_ani;
extern std::vector<double> Bx_exc, By_exc, Bz_exc;
extern std::vector<double> Bx_lon, By_lon, Bz_lon;
extern std::vector<double> Bx_eff, By_eff, Bz_eff;

extern std::vector<double> torque_x, torque_y, torque_z, torque_mod;

//---Functions
int uniax_anis_f(int n_cells,
				 std::vector<double>K, std::vector<double>Ms,
				 std::vector<double>ex, std::vector<double> ey, std::vector<double> ez,
				 std::vector<double>mx, std::vector<double> my, std::vector<double> mz,
				 std::vector<double>&Bx_ani, std::vector<double> &By_ani, std::vector<double> &Bz_ani,
				 std::vector<int>mat_id);
				

int zeeman_f(int n_cells,
			 double B_app, double bx, double by, double bz,
			 std::vector<double> &Bx_app, std::vector<double> &By_app, std::vector<double> &Bz_app);

int exchange_f(int n_cells, double lengthscale,
			   std::vector<double> m_e, std::vector<double> Ms0_SI, std::vector<unsigned long long int> macrocell_size, std::vector<std::vector<double>>A_T_matrix,
			   std::vector<int> int_list, std::vector<int> start_neighbours, std::vector<int> end_neighbours,
			   std::vector<int> material_id,
			   std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
			   std::vector<double> &Bx_exc, std::vector<double> &By_exc, std::vector<double> &Bz_exc);


int longitudinal_f(int n_cells,
				   std::vector<double> chi_par, std::vector<double> m_e,
				   std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
				   std::vector<double> &Bx_lon, std::vector<double> &By_lon, std::vector<double> &Bz_lon,
				   std::vector<int>mat_id);

int effective_f(int n_cells,
				std::vector<double> Bx_ani, std::vector<double> By_ani, std::vector<double> Bz_ani,
				std::vector<double> Bx_app, std::vector<double> By_app, std::vector<double> Bz_app,
				std::vector<double> Bx_exc, std::vector<double> By_exc, std::vector<double> Bz_exc,
				std::vector<double> Bx_lon, std::vector<double> By_lon, std::vector<double> Bz_lon,
				std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff);

int effective_torque_f(int n_cells,
					   std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
					   std::vector<double> Bx_eff, std::vector<double> By_eff, std::vector<double> Bz_eff,
					   std::vector<double> &torque_x, std::vector<double> &torque_y, std::vector<double> &torque_z,
					   std::vector<double>  &torque_mod);
int calculate();



}

namespace field{

	namespace internal{

		int alloc_memory(int n_cells);


	}

}

#endif
//---End of field.h file.