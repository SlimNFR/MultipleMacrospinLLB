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
    int id; //stores the material id of the cell
    double uniax_energy = 0.0; // this will be returned at the end

    for(int i=0; i<n_cells; i++)
    {//loop cells
        id=mat_id[i];
        uniax_energy -= K[id]*V[id]*pow((mx[i]*ex[id] + my[i]*ey[id] + mz[i]*ez[id]),2.0);

    }
    
    return uniax_energy; //Uni-axial magnetocrystalline anisotropy energy: [J]
 }


 double zeeman_f(int n_cells,
                 double B, double bx, double by, double bz,
                 std::vector<double> Ms, std::vector<double> V,
                 std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                 std::vector<double> mat_id)
 {
   int id; //stores the material id of the cell
   double zeeman_energy = 0.0;//this will be returned at the end

   for(int i=0; i<n_cells; i++)
   {

      id=mat_id[i];
      zeeman_energy-= V[id]*Ms[id]*B*(mx[i]*bx + my[i]*by + mz[i]*bz);

   }

 	return zeeman_energy; //Zeeman energy: [J]

 }


}


//---End of energy.cpp file.

