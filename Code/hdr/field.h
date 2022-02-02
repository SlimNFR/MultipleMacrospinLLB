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

extern std::vector<double> torque_app_x, torque_app_y, torque_app_z, torque_app_mod;
extern std::vector<double> torque_ani_x, torque_ani_y, torque_ani_z, torque_ani_mod;
extern std::vector<double> torque_exc_x, torque_exc_y, torque_exc_z, torque_exc_mod;
extern std::vector<double> torque_lon_x, torque_lon_y, torque_lon_z, torque_lon_mod;
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
			   std::vector<double> m_e, std::vector<double> Ms0_SI, std::vector<unsigned long long int> macrocell_size,
			   std::vector<unsigned long long int>unitcell_size, std::vector<std::vector<double>>A_T_matrix,
			   std::vector<int> int_list, std::vector<int> start_neighbours, std::vector<int> end_neighbours,
			   std::vector<int> material_id,
			   std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
			   std::vector<double> &Bx_exc, std::vector<double> &By_exc, std::vector<double> &Bz_exc);


int longitudinal_f(int n_cells,
				   double T,std::vector<double>Tc,
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

int torque_app_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_app, std::vector<double> By_app, std::vector<double> Bz_app,
				 std::vector<double> &torque_app_x, std::vector<double> &torque_app_y, std::vector<double> &torque_app_z,
			     std::vector<double> &torque_app_mod);

int torque_ani_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_ani, std::vector<double> By_ani, std::vector<double> Bz_ani,
				 std::vector<double> &torque_ani_x, std::vector<double> &torque_ani_y, std::vector<double> &torque_ani_z,
			     std::vector<double> &torque_ani_mod);

int torque_exc_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_exc, std::vector<double> By_exc, std::vector<double> Bz_exc,
				 std::vector<double> &torque_exc_x, std::vector<double> &torque_exc_y, std::vector<double> &torque_exc_z,
			     std::vector<double> &torque_exc_mod);

int torque_lon_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_lon, std::vector<double> By_lon, std::vector<double> Bz_lon,
				 std::vector<double> &torque_lon_x, std::vector<double> &torque_lon_y, std::vector<double> &torque_lon_z,
			     std::vector<double> &torque_lon_mod);

int effective_torque_f(int n_cells,
					   bool remove_precession_term,
					   bool remove_longitudin_term,
					   bool remove_transverse_term,
					   std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
					   std::vector<double> Bx_eff, std::vector<double> By_eff, std::vector<double> Bz_eff,
					   std::vector<double> &torque_x, std::vector<double> &torque_y, std::vector<double> &torque_z,
					   std::vector<double>  &torque_mod);

int adjust_field_f(int n_cells,
				   bool force_DW_formation,
				   std::vector<int>left_edge_spins,
				   std::vector<int>right_edge_spins,
				   std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff);

int adjust_torque_f(int n_cells,
					bool force_DW_formation,
					std::vector<int>left_edge_spins,
					std::vector<int>right_edge_spins,
					std::vector<double> &torque_x, std::vector<double> &torque_y, std::vector<double> &torque_z,
					std::vector<double> &torque_mod);


int calculate();



}

namespace field{

	namespace internal{

		int alloc_memory(int n_cells);


	}

}

#endif
//---End of field.h file.