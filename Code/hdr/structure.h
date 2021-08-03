//---This is the structure.h file.
#pragma once

#ifndef STRUCTURE_H
#define STRUCTURE_H

//---Standard libraries
#include<vector>

//---User-defined libraries
#include"input.h"

//---Namespace macrospin

namespace macrospin{
    
//---Variables
extern std::vector<double>mx_0,my_0,mz_0;
extern std::vector<double>mx,my,mz;


//---Functions
int update_magnetisation(double mx_in,double my_in, double mz_in,
						 double &mx, double &my, double &mz);

}

namespace material{

//---Variables

extern std::vector<double> xcoord;
extern std::vector<double> ycoord;
extern std::vector<double> zcoord;
extern std::vector<int> id;

extern std::vector<int> interaction_list;
extern std::vector<int> start_neighbours;
extern std::vector<int> end_neighbours;



//---Functions

int create_interaction_list(int n_cells,
                            std::vector<double> xcoord,
                            std::vector<double> ycoord,
                            std::vector<double> zcoord,
                            std::vector<int> &int_list,
                            std::vector<int> &start,
                            std::vector<int> &end);


int generate_crystal_structure_f(int n_materials,
                                 int n_cells,
                                 std::vector<int>nx, std::vector<int>ny, std::vector<int>nz,
                                 std::vector<int> &mat_id,
                                 std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);




}


namespace macrospin{


    namespace  internal{


        int alloc_memory(int n_cells);
        int set_initial_config(int n_cells,
                           std::vector<double> mx0_in,std::vector<double> my0_in,std::vector<double> mz0_in,
                           std::vector<double> &mx0_out,std::vector<double> &my0_out,std::vector<double> &mz0_out,
                           std::vector<int>material_id);
    }
}

#endif
//---End of structure.h file.