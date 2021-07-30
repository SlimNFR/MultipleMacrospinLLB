//---This is the energy.cpp file. It defines the energy.h functions.


//---Standard libraries
#include<cmath>
#include<vector>

//---User-defined libraries
#include"energy.h"
#include"input.h"

//---Namespace energy

namespace energy{
	
//---Variables
double uniax_anis=0.0;
double zeeman=0.0;


//---Functions 
 double uniax_anis_f(int n_cells,
                     std::vector<double> K, std::vector<double> V,
                     std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                     std::vector<double> ex, std::vector<double> ey, std::vector<double> ez,
                     std::vector<int> mat_id)
 {
    int mat; //stores the material id of the cell
    double uniax_energy = 0.0; // this will be returned at the end

    for(int cell=0; cell<n_cells; cell++)
    {//loop cells
        mat=mat_id[cell];
        uniax_energy -= K[mat]*V[mat]*pow((mx[cell]*ex[mat] + my[cell]*ey[mat] + mz[cell]*ez[mat]),2.0); //E_ani : [J]

    }
    
    return uniax_energy; //Uni-axial magnetocrystalline anisotropy energy: [J]
 }


 double zeeman_f(int n_cells,
                 double B, double bx, double by, double bz,
                 std::vector<double> Ms, std::vector<double> V,
                 std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                 std::vector<double> mat_id)
 {
   int mat; //stores the material id of the cell
   double zeeman_energy = 0.0;//this will be returned at the end

   for(int cell=0; cell<n_cells; cell++)
   {

      mat=mat_id[cell];
      zeeman_energy-= V[mat]*Ms[mat]*B*(mx[cell]*bx + my[cell]*by + mz[cell]*bz); //E_zeeman : [J]

   }

 	return zeeman_energy; //Zeeman energy: [J]

 }


}


//---End of energy.cpp file.

