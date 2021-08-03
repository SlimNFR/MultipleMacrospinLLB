//---This is the field.cpp file. It defines the field.h functions.


//---Standard libraries
#include<iostream>
#include<cmath>

//---User-defined libraries
#include"field.h"
#include"structure.h"
#include"input.h"

//---Namespace field

namespace field{
	
//---Variables
std::vector<double> Bx_app, By_app, Bz_app;
std::vector<double> Bx_ani, By_ani, Bz_ani;
std::vector<double> Bx_exc, By_exc, Bz_exc;
std::vector<double> Bx_lon, By_lon, Bz_lon;
std::vector<double> Bx_eff, By_eff, Bz_eff;

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
			   std::vector<double> m_e, std::vector<double> Ms0_SI, std::vector<int> macrocell_size, std::vector<std::vector<double>>A_T_matrix,
			   std::vector<int> int_list, std::vector<int> start_neighbours, std::vector<int> end_neighbours,
			   std::vector<int> material_id,
			   std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
			   std::vector<double> &Bx_exc, std::vector<double> &By_exc, std::vector<double> &Bz_exc)
{

	for(int cell=0; cell<n_cells; cell++) //Loop each cell
	{
		int mat_id_cell=material_id[cell]; //Get material id of cell
		double pre_factor = (pow(m_e[mat_id_cell],2.0)*Ms0_SI[mat_id_cell]*pow(macrocell_size[mat_id_cell]*lengthscale, 2.0));
		for(int count_neighbour=start_neighbours[cell]; count_neighbour<end_neighbours[cell]; count_neighbour++) //Count all cell's neighbours
		{
			int neighbour = int_list[count_neighbour]; //Get neighbour from interaction list
			int mat_id_neighbour=material_id[neighbour]; //Get material id of neighbour
			double A = A_T_matrix[mat_id_cell][mat_id_neighbour]; //Get exchange from exchange matrix
		//	std::cout<<"A: "<<A<<" |prefactor: "<<pre_factor<<"\n";
			Bx_exc[cell] -= (2.0*A/pre_factor)*(mx[neighbour] - mx[cell]); //  [T]
			By_exc[cell] -= (2.0*A/pre_factor)*(my[neighbour] - my[cell]);
			Bz_exc[cell] -= (2.0*A/pre_factor)*(mz[neighbour] - mz[cell]);

		}

		//std::cout<<"Cell: "<<cell<<" |Bx_exc_cell:"<<Bx_exc[cell]<<" |By__exc_cell:"<<By_exc[cell]<<" |Bz__exc_cell:"<<Bz_exc[cell]<<"\n";

	}


	return 0;
}


int longitudinal_f(int n_cells,
				   std::vector<double> chi_par, std::vector<double> m_e,
				   std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
				   std::vector<double> &Bx_lon, std::vector<double> &By_lon, std::vector<double> &Bz_lon,
				   std::vector<int>mat_id)
{
	//calculates the longitudinal field components
	//I assume I work below Tc for now
	double m_squared; //This will save the value of mxmx + mymy + mzmz
	double pre_factor; //This is the prefactor in front of the long. field (1/(2chi_par))(1-m^2/me^2)
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
		pre_factor = (0.5*(1.0/chi_par[mat])) * (1.0 - (m_squared/(m_e[mat]*m_e[mat]))); 
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

	for(int cell=0; cell<n_cells;cell++)
	{
		Bx_eff[cell] = Bx_ani[cell] + Bx_app[cell] + Bx_exc[cell] + Bx_lon[cell]; //B_eff: [T]
		By_eff[cell] = By_ani[cell] + By_app[cell] + By_exc[cell] + By_lon[cell];
		Bz_eff[cell] = Bz_ani[cell] + Bz_app[cell] + Bz_exc[cell] + Bz_lon[cell];
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

	torque_x[cell] = (my[cell]*Bz_eff[cell] - mz[cell]*By_eff[cell]);
	torque_y[cell] = (mz[cell]*Bx_eff[cell] - mx[cell]*Bz_eff[cell]);
	torque_z[cell] = (mx[cell]*By_eff[cell] - my[cell]*Bx_eff[cell]);

	torque_mod[cell] = sqrt(torque_x[cell]*torque_x[cell]+
				 	  	    torque_y[cell]*torque_y[cell]+
				 	  		torque_z[cell]*torque_z[cell]);	

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
					  input::m_e, input::Ms0_SI, input::unitcell_size, input::A_T_matrix,
					  material::interaction_list, material::start_neighbours, material::end_neighbours, material::id, 
					  macrospin::mx, macrospin::my, macrospin::mz,
					  field::Bx_exc, field::By_exc, field::Bz_exc);

	field::longitudinal_f(input::n_cells,
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

	field::effective_torque_f(input::n_cells,
							  macrospin::mx, macrospin::my, macrospin::mz,
							  field::Bx_eff, field::By_eff, field::Bz_eff,
							  field::torque_x, field::torque_y, field::torque_z,
							  field::torque_mod);

	

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
			

			return 0;
		}
	}
}


//---End of field.cpp file.

