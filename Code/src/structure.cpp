//---This is the structure.cpp file. It defines the structure.h functions.


//---Standard libraries
#include<vector>
#include<iostream>
//---User-defined libraries
#include"structure.h"
#include"input.h"
//---Namespace macrospin

namespace macrospin{
	
//---Variables
std::vector<double>mx(input::n_cells);
std::vector<double>my(input::n_cells);
std::vector<double>mz(input::n_cells);

//---Functions
int update_magnetisation(double mx_in,double my_in, double mz_in,
			 double &mx, double &my, double &mz)
{	//self-explanatory
	mx = mx_in; my = my_in; mz = mz_in;
	return 0;
}


}

//---Namespace material
namespace material{

//---Variables

//These four arrays will save the real crystallographic positions of the macrospins and their material id.
std::vector<double> xcoord; 
std::vector<double> ycoord;
std::vector<double> zcoord;
std::vector<int> id;

//---Functions

void generate_f(int n_materials,
                int n_cells,
                std::vector<int>nx, std::vector<int>ny, std::vector<int>nz,
                std::vector<int> &mat_id,
                std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
{//This function will assign real coordinates and a material id to each of the macrospins in my system.

  //I'll assume in terms of widths and length the materials are the same!
  
  mat_id.resize(n_cells);
  x.resize(n_cells);
  y.resize(n_cells);
  z.resize(n_cells);

  int count_cell = 0;
  int count_material = 0;
  int base=0;
  int mat = 0;

  //---set up the real coordinates of each macrospin/cell
  //---assign each cell a material id

  for(mat = 0; mat<n_materials;mat++)
  {
    for(int idz=0; idz<nz[mat]; idz++) 
    {
     for(int idy=0; idy<ny[mat]; idy++)
     {
      for(int idx=0; idx<nx[mat]; idx++)
      {
       
       //for current mat, start counting the atoms from the (0,0,height_previous_material)
       x[count_cell] = idx;
       y[count_cell] = idy;
       z[count_cell] = idz+base; 
       
       mat_id[count_cell] = mat;
       std::cout<<x[count_cell]<<" "<<y[count_cell]<<" "<<z[count_cell]<<" "<<mat_id[count_cell]<<"\n";

       
       count_cell++;
      
      }
     }
    }
    base += nz[mat];
  }  
}
}

//---End of structure.cpp file.

