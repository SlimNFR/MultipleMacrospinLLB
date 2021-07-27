//---This is the main.cpp file. It calls the main routines of the code.


//---Standard libraries
#include<iostream>

//---User-defined libraries
#include"input.h"
#include"structure.h"

//#include"solver.h"
//#include"field.h"
//#include"init.h"
//#include"tempscaling.h"
//#include"interpolate.h"
//#include"output.h"
//#include"simulation.h"

int main()
{

input::read_material_parameters();
input::read_simulation_parameters();
material::generate_f(input::n_materials,
					 input::n_cells,
					 input::nx_cells, input::ny_cells, input::nz_cells,
					 material::id,
					 material::xcoord, material::ycoord, material::zcoord);

//init::internal::run();

//init::internal::end();

return 0;
}





//---End of main.cpp file
