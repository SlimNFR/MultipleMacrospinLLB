//---This is the field.cpp file. It defines the field.h functions.


//---Standard libraries
#include<iostream>
#include<cmath>

//---User-defined libraries
#include"field.h"
#include"structure.h"
#include"input.h"
#include"utils.h"

//---Namespace field

namespace field{
	
//---Variables
std::vector<double> Bx_app, By_app, Bz_app;
std::vector<double> Bx_ani, By_ani, Bz_ani;
std::vector<double> Bx_exc, By_exc, Bz_exc;
std::vector<double> Bx_lon, By_lon, Bz_lon;
std::vector<double> Bx_eff, By_eff, Bz_eff;

std::vector<double> torque_app_x, torque_app_y, torque_app_z, torque_app_mod;
std::vector<double> torque_ani_x, torque_ani_y, torque_ani_z, torque_ani_mod;
std::vector<double> torque_exc_x, torque_exc_y, torque_exc_z, torque_exc_mod;
std::vector<double> torque_lon_x, torque_lon_y, torque_lon_z, torque_lon_mod;
std::vector<double> torque_x, torque_y, torque_z, torque_mod;

//---Functions
int uniax_anis_f(int n_cells,
				 std::vector<double>K, std::vector<double>Ms,
				 std::vector<double>ex, std::vector<double> ey, std::vector<double> ez,
				 std::vector<double>mx, std::vector<double> my, std::vector<double> mz,
				 std::vector<double>&Bx_ani, std::vector<double> &By_ani, std::vector<double> &Bz_ani,
				 std::vector<int>mat_id)
{
	//calculates the anisotropy field components
	int mat; //will store the material id.

	for(int cell=0; cell<n_cells; cell++)
	{	
		mat = mat_id[cell];//Get the material id of the spin
		double Hk = 2.0*K[mat]/Ms[mat];
		Bx_ani[cell] = Hk*(mx[cell]*ex[mat] + my[cell]*ey[mat] + mz[cell]*ez[mat])*ex[mat];	//B_ani : [T]
		By_ani[cell] = Hk*(mx[cell]*ex[mat] + my[cell]*ey[mat] + mz[cell]*ez[mat])*ey[mat];
		Bz_ani[cell] = Hk*(mx[cell]*ex[mat] + my[cell]*ey[mat] + mz[cell]*ez[mat])*ez[mat];
		//std::cout<<"Anisotropy is: "<<K[mat]<<"\n";
		//std::cout<<"Hk: "<<Hk<<" Ms:"<<Ms[mat]<<"K: "<<K[mat]<<"\n";
		//std::cout<<"Bx_ani: "<<Bx_ani[cell]<<"By_ani: "<<By_ani[cell]<<"Bz_ani: "<<Bz_ani[cell]<<"\n";
	}

	return 0;
}
				

int zeeman_f(int n_cells,
			 double B_app, double bx, double by, double bz,
			 std::vector<double> &Bx_app, std::vector<double> &By_app, std::vector<double> &Bz_app)
{
	//calculates the zeeman field components
	
	for(int cell=0; cell<n_cells; cell++)
	{
		Bx_app[cell]=B_app*bx; //B_app: [T]
		By_app[cell]=B_app*by;
		Bz_app[cell]=B_app*bz;

		//std::cout<<"Bx_app: "<<Bx_app[cell]<<"By_app: "<<By_app[cell]<<"Bz_app: "<<Bz_app[cell]<<"\n";

	}
	return 0;
}



int exchange_f(int n_cells, double lengthscale,
			   std::vector<double> m_e, std::vector<double> Ms0_SI, std::vector<unsigned long long int> macrocell_size, std::vector<std::vector<double>>A_T_matrix,
			   std::vector<int> int_list, std::vector<int> start_neighbours, std::vector<int> end_neighbours,
			   std::vector<int> material_id,
			   std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
			   std::vector<double> &Bx_exc, std::vector<double> &By_exc, std::vector<double> &Bz_exc)
{
	//I need to reset the exchange field everytime.
	std::fill(Bx_exc.begin(), Bx_exc.end(), 0.0);
	std::fill(By_exc.begin(), By_exc.end(), 0.0);
	std::fill(Bz_exc.begin(), Bz_exc.end(), 0.0);

	int neighbour, mat_id_neighbour, mat_id_cell;
	double A;
	for(int cell=0; cell<n_cells; cell++) //Loop each cell
	{
		mat_id_cell=material_id[cell]; //Get material id of cell
		double pre_factor = (pow(m_e[mat_id_cell],2.0)*Ms0_SI[mat_id_cell]*pow(macrocell_size[mat_id_cell]*lengthscale, 2.0));
		for(int count_neighbour=start_neighbours[cell]; count_neighbour<end_neighbours[cell]; count_neighbour++) //Count all cell's neighbours
		{

			neighbour = int_list[count_neighbour]; //Get neighbour from interaction list
			mat_id_neighbour=material_id[neighbour]; //Get material id of neighbour
			A = A_T_matrix[mat_id_cell][mat_id_neighbour]; //Get exchange from exchange matrix
			//std::cout<<"Exchange is: "<<A<<"\n";
			//Calculate exchange
			Bx_exc[cell] += (2.0*A/pre_factor)*(mx[neighbour] - mx[cell]); //  [T]
			By_exc[cell] += (2.0*A/pre_factor)*(my[neighbour] - my[cell]);
			Bz_exc[cell] += (2.0*A/pre_factor)*(mz[neighbour] - mz[cell]);
		}

		//std::cout<<"Cell: "<<cell<<" |Bx_exc_cell:"<<Bx_exc[cell]<<" |By__exc_cell:"<<By_exc[cell]<<" |Bz__exc_cell:"<<Bz_exc[cell]<<"\n";
				 
	}


	return 0;
}


int longitudinal_f(int n_cells,
				   double T,std::vector<double>Tc,
				   std::vector<double> chi_par, std::vector<double> m_e,
				   std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
				   std::vector<double> &Bx_lon, std::vector<double> &By_lon, std::vector<double> &Bz_lon,
				   std::vector<int>mat_id)
{
	//calculates the longitudinal field components

	double m_squared; //This will save the value of mxmx + mymy + mzmz
	double pre_factor; //This is the prefactor in front of the long. field 
	int mat; //This will store the material id
	if(input::T==0)
	{	//If the temperature is zero then no longitudinal field should be present..
		std::fill(Bx_lon.begin(), Bx_lon.end(), 0.0);
		std::fill(By_lon.begin(), By_lon.end(), 0.0);
		std::fill(Bz_lon.begin(), Bz_lon.end(), 0.0);
		return 0;
	}
	
	for(int cell=0; cell<n_cells; cell++)
	{
		mat = mat_id[cell];//Get material id.
		m_squared = mx[cell]*mx[cell] + my[cell]*my[cell] + mz[cell]*mz[cell]; 
		if(T<Tc[mat])pre_factor = (0.5*(1.0/chi_par[mat])) * (1.0 - (m_squared/(m_e[mat]*m_e[mat]))); 
		else if(T>Tc[mat])pre_factor= -1.0/chi_par[mat]*(1.0 - (3.0*Tc[mat]/(5.0*(T-Tc[mat])))*m_squared);
		Bx_lon[cell] = pre_factor*mx[cell]; //B_lon: [T]
		By_lon[cell] = pre_factor*my[cell];
		Bz_lon[cell] = pre_factor*mz[cell];

		//std::cout<<"Bx_lon: "<<Bx_lon[cell]<<"By_lon: "<<By_lon[cell]<<"Bz_lon: "<<Bz_lon[cell]<<"\n";
	}

	return 0;

}




int effective_f(int n_cells,
				std::vector<double> Bx_ani, std::vector<double> By_ani, std::vector<double> Bz_ani,
				std::vector<double> Bx_app, std::vector<double> By_app, std::vector<double> Bz_app,
				std::vector<double> Bx_exc, std::vector<double> By_exc, std::vector<double> Bz_exc,
				std::vector<double> Bx_lon, std::vector<double> By_lon, std::vector<double> Bz_lon,
				std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff)
{
	//calculates the total field (effective) components

	//std::vector<double>B_mod;
	//B_mod.resize(n_cells,0.0);
	//double max_Bmod;
	//int id_max_Bmod;

	for(int cell=0; cell<n_cells;cell++)
	{
		Bx_eff[cell] = Bx_ani[cell] + Bx_app[cell] + Bx_exc[cell] + Bx_lon[cell]; //B_eff: [T]
		By_eff[cell] = By_ani[cell] + By_app[cell] + By_exc[cell] + By_lon[cell];
		Bz_eff[cell] = Bz_ani[cell] + Bz_app[cell] + Bz_exc[cell] + Bz_lon[cell];

		/*
		B_mod[cell] = sqrt(Bx_eff[cell]*Bx_eff[cell] +
						   By_eff[cell]*By_eff[cell] +
						   Bz_eff[cell]*Bz_eff[cell]);
		*/

	}

	//utils::max_element_1D_vec(B_mod, max_Bmod, id_max_Bmod, false);
	//std::cout<<"Cell max B_mod "<<id_max_Bmod<<" |max Bmod: "<<max_Bmod<<"\n";


	return 0;

}


int torque_app_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_app, std::vector<double> By_app, std::vector<double> Bz_app,
				 std::vector<double> &torque_app_x, std::vector<double> &torque_app_y, std::vector<double> &torque_app_z,
			     std::vector<double> &torque_app_mod)
{

	for(int cell=0; cell<n_cells; cell++)
	{



	torque_app_x[cell] = (my[cell]*Bz_app[cell] - mz[cell]*By_app[cell]);
	torque_app_y[cell] = (mz[cell]*Bx_app[cell] - mx[cell]*Bz_app[cell]);
	torque_app_z[cell] = (mx[cell]*By_app[cell] - my[cell]*Bx_app[cell]);

	torque_app_mod[cell] = sqrt(torque_app_x[cell]*torque_app_x[cell]+
				 	  	    	torque_app_y[cell]*torque_app_y[cell]+
				 	  			torque_app_z[cell]*torque_app_z[cell]);	

	}

	return 0;
}


int torque_ani_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_ani, std::vector<double> By_ani, std::vector<double> Bz_ani,
				 std::vector<double> &torque_ani_x, std::vector<double> &torque_ani_y, std::vector<double> &torque_ani_z,
			     std::vector<double> &torque_ani_mod)
{

	for(int cell=0; cell<n_cells; cell++)
	{


	torque_ani_x[cell] = (my[cell]*Bz_ani[cell] - mz[cell]*By_ani[cell]);
	torque_ani_y[cell] = (mz[cell]*Bx_ani[cell] - mx[cell]*Bz_ani[cell]);
	torque_ani_z[cell] = (mx[cell]*By_ani[cell] - my[cell]*Bx_ani[cell]);

	torque_ani_mod[cell] = sqrt(torque_ani_x[cell]*torque_ani_x[cell]+
				 	  	    	torque_ani_y[cell]*torque_ani_y[cell]+
				 	  			torque_ani_z[cell]*torque_ani_z[cell]);	

	}

	return 0;
}

int torque_exc_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_exc, std::vector<double> By_exc, std::vector<double> Bz_exc,
				 std::vector<double> &torque_exc_x, std::vector<double> &torque_exc_y, std::vector<double> &torque_exc_z,
			     std::vector<double> &torque_exc_mod)
{

	for(int cell=0; cell<n_cells; cell++)
	{


	torque_exc_x[cell] = (my[cell]*Bz_exc[cell] - mz[cell]*By_exc[cell]);
	torque_exc_y[cell] = (mz[cell]*Bx_exc[cell] - mx[cell]*Bz_exc[cell]);
	torque_exc_z[cell] = (mx[cell]*By_exc[cell] - my[cell]*Bx_exc[cell]);

	torque_exc_mod[cell] = sqrt(torque_exc_x[cell]*torque_exc_x[cell]+
				 	  	    	torque_exc_y[cell]*torque_exc_y[cell]+
				 	  			torque_exc_z[cell]*torque_exc_z[cell]);	

	}

	return 0;
}

int torque_lon_f(int n_cells,
				 std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
				 std::vector<double> Bx_lon, std::vector<double> By_lon, std::vector<double> Bz_lon,
				 std::vector<double> &torque_lon_x, std::vector<double> &torque_lon_y, std::vector<double> &torque_lon_z,
			     std::vector<double> &torque_lon_mod)
{

	for(int cell=0; cell<n_cells; cell++)
	{



	torque_lon_x[cell] = (my[cell]*Bz_lon[cell] - mz[cell]*By_lon[cell]);
	torque_lon_y[cell] = (mz[cell]*Bx_lon[cell] - mx[cell]*Bz_lon[cell]);
	torque_lon_z[cell] = (mx[cell]*By_lon[cell] - my[cell]*Bx_lon[cell]);

	torque_lon_mod[cell] = sqrt(torque_lon_x[cell]*torque_lon_x[cell]+
				 	  	    	torque_lon_y[cell]*torque_lon_y[cell]+
				 	  			torque_lon_z[cell]*torque_lon_z[cell]);	

	}

	return 0;
}




int effective_torque_f(int n_cells,
					   std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
					   std::vector<double> Bx_eff, std::vector<double> By_eff, std::vector<double> Bz_eff,
					   std::vector<double> &torque_x, std::vector<double> &torque_y, std::vector<double> &torque_z,
					   std::vector<double> &torque_mod)
{

	for(int cell=0; cell<n_cells; cell++)
	{


	double m_mod, B_eff_mod;
	m_mod=sqrt(mx[cell]*mx[cell] + my[cell]*my[cell] + mz[cell]*mz[cell]);


	/*For NORMAL Version mxH comment the 4 next lines*/

	Bx_eff[cell]=Bx_eff[cell]-field::Bx_lon[cell];
	By_eff[cell]=By_eff[cell]-field::By_lon[cell];
	Bz_eff[cell]=Bz_eff[cell]-field::Bz_lon[cell];
	B_eff_mod=sqrt(Bx_eff[cell]*Bx_eff[cell] + By_eff[cell]*By_eff[cell] + Bz_eff[cell]*Bz_eff[cell]);

	
	
	torque_x[cell] = (my[cell]*Bz_eff[cell] - mz[cell]*By_eff[cell])/(B_eff_mod*m_mod);
	torque_y[cell] = (mz[cell]*Bx_eff[cell] - mx[cell]*Bz_eff[cell])/(B_eff_mod*m_mod);
	torque_z[cell] = (mx[cell]*By_eff[cell] - my[cell]*Bx_eff[cell])/(B_eff_mod*m_mod);
	
	
	torque_mod[cell] = sqrt(torque_x[cell]*torque_x[cell]+
				 	  	    torque_y[cell]*torque_y[cell]+
				 	  		torque_z[cell]*torque_z[cell]);	
	
	}

	return 0;
}

int adjust_field_f(int n_cells,
				   bool force_DW_formation,
				   std::vector<double> &Bx_eff, std::vector<double> &By_eff, std::vector<double> &Bz_eff)
{	//this function helps me readjust torques 

	if(force_DW_formation==true)
	{	//in the case of DW formation I want the effective fields on the edge to be zero.
		Bx_eff[0]=By_eff[0]=Bz_eff[0]=0.0;
		Bx_eff[n_cells-1]=By_eff[n_cells-1]=Bz_eff[n_cells-1]=0.0;
	}

	return 0;
}

int adjust_torque_f(int n_cells,
					bool force_DW_formation,
					std::vector<double> &torque_x, std::vector<double> &torque_y, std::vector<double> &torque_z,
					std::vector<double> &torque_mod)
{	//this function helps me readjust torques 

	if(force_DW_formation==true)
	{	//in the case of DW formation I want the torques on the edge spins to be zero.
		torque_mod[0]=0.0;
		torque_mod[n_cells-1]=0.0;
	}

	return 0;
}

int calculate()
{
	field::uniax_anis_f(input::n_cells,
						input::K_T, input::Ms_T,
						input::ex, input::ey, input::ez,
						macrospin::mx, macrospin::my, macrospin::mz,
						field::Bx_ani, field::By_ani, field::Bz_ani,
						material::id);

	field::zeeman_f(input::n_cells,
					input::B_app, input::bx, input::by, input::bz,
					field::Bx_app, field::By_app, field::Bz_app);

	
	field::exchange_f(input::n_cells,input::lengthscale,
					  input::m_e, input::Ms0_SI, input::macrocell_size, input::A_T_matrix,
					  material::interaction_list, material::start_neighbours, material::end_neighbours, material::id, 
					  macrospin::mx, macrospin::my, macrospin::mz,
					  field::Bx_exc, field::By_exc, field::Bz_exc);

	field::longitudinal_f(input::n_cells,
						  input::T, input::Tc,
						  input::chi_par, input::m_e,
						  macrospin::mx, macrospin::my, macrospin::mz,
						  field::Bx_lon, field::By_lon, field::Bz_lon,
						  material::id);

	field::effective_f(input::n_cells,
					   field::Bx_ani, field::By_ani, field::Bz_ani,
					   field::Bx_app, field::By_app, field::Bz_app,
					   field::Bx_exc, field::By_exc, field::Bz_exc,
					   field::Bx_lon, field::By_lon, field::Bz_lon,
					   field::Bx_eff, field::By_eff, field::Bz_eff);
	
	field::adjust_field_f(input::n_cells,
						  input::force_DW_formation,
						  field::Bx_eff, field::By_eff, field::Bz_eff);
	
	field::torque_app_f(input::n_cells,
						macrospin::mx, macrospin::my, macrospin::mz,
						field::Bx_app, field::By_app, field::Bz_app,
						field::torque_app_x, field::torque_app_y, field::torque_app_z,
						field::torque_app_mod);


	field::torque_ani_f(input::n_cells,
						macrospin::mx, macrospin::my, macrospin::mz,
						field::Bx_ani, field::By_ani, field::Bz_ani,
						field::torque_ani_x, field::torque_ani_y, field::torque_ani_z,
						field::torque_ani_mod);


	field::torque_exc_f(input::n_cells,
						macrospin::mx, macrospin::my, macrospin::mz,
						field::Bx_exc, field::By_exc, field::Bz_exc,
						field::torque_exc_x, field::torque_exc_y, field::torque_exc_z,
						field::torque_exc_mod);

	field::torque_lon_f(input::n_cells,
						macrospin::mx, macrospin::my, macrospin::mz,
						field::Bx_lon, field::By_lon, field::Bz_lon,
						field::torque_lon_x, field::torque_lon_y, field::torque_lon_z,
						field::torque_lon_mod);

	field::effective_torque_f(input::n_cells,
							  macrospin::mx, macrospin::my, macrospin::mz,
							  field::Bx_eff, field::By_eff, field::Bz_eff,
							  field::torque_x, field::torque_y, field::torque_z,
							  field::torque_mod);
	

	

	/*
	field::adjust_torque_f(input::n_cells,
						   input::force_DW_formation,
						   field::torque_x, field::torque_y, field::torque_z,
						   field::torque_mod);
	*/

	//adjust::field OR 
	//adjust::torque


	

	return 0;
}


}

namespace field{


	namespace internal{


		int alloc_memory(int n_cells)
		{ //This function allocates memory for the

			field::Bx_app.resize(n_cells,0.0);
			field::By_app.resize(n_cells,0.0);
			field::Bz_app.resize(n_cells,0.0);
			field::Bx_ani.resize(n_cells,0.0);
			field::By_ani.resize(n_cells,0.0);
			field::Bz_ani.resize(n_cells,0.0);
			field::Bx_exc.resize(n_cells,0.0);
			field::By_exc.resize(n_cells,0.0);
			field::Bz_exc.resize(n_cells,0.0);
			field::Bx_lon.resize(n_cells,0.0);
			field::By_lon.resize(n_cells,0.0);
			field::Bz_lon.resize(n_cells,0.0);
			field::Bx_eff.resize(n_cells,0.0);
			field::By_eff.resize(n_cells,0.0);
			field::Bz_eff.resize(n_cells,0.0);
			field::torque_x.resize(n_cells,0.0);
			field::torque_y.resize(n_cells,0.0);
			field::torque_z.resize(n_cells,0.0);
			field::torque_mod.resize(n_cells,0.0);
			field::torque_app_x.resize(n_cells,0.0);
			field::torque_app_y.resize(n_cells,0.0);
			field::torque_app_z.resize(n_cells,0.0);
			field::torque_app_mod.resize(n_cells,0.0);
			field::torque_ani_x.resize(n_cells,0.0);
			field::torque_ani_y.resize(n_cells,0.0);
			field::torque_ani_z.resize(n_cells,0.0);
			field::torque_ani_mod.resize(n_cells,0.0);
			field::torque_exc_x.resize(n_cells,0.0);
			field::torque_exc_y.resize(n_cells,0.0);
			field::torque_exc_z.resize(n_cells,0.0);
			field::torque_exc_mod.resize(n_cells,0.0);
			field::torque_lon_x.resize(n_cells,0.0);
			field::torque_lon_y.resize(n_cells,0.0);
			field::torque_lon_z.resize(n_cells,0.0);
			field::torque_lon_mod.resize(n_cells,0.0);
			

			return 0;
		}
	}
}


//---End of field.cpp file.

