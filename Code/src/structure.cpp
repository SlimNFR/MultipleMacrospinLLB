//---This is the structure.cpp file. It defines the structure.h functions.


//---Standard libraries
#include<vector>
#include<iostream>
#include<deque>
#include<cmath>
#include<fstream>
#include<iostream>
#include<string>
//---User-defined libraries
#include"structure.h"

//---Namespace macrospin

namespace macrospin{
	
//---Variables
std::vector<double>mx_0,my_0,mz_0;
std::vector<double>mx,my,mz;


//---Functions
int update_magnetisation(double mx_in,double my_in, double mz_in,
			 double &mx, double &my, double &mz)
{	//self-explanatory
	mx = mx_in; my = my_in; mz = mz_in;
	return 0;
}

int read_config_from_file(int n_cells,
                          std::vector<double> &mx0_out,std::vector<double> &my0_out,std::vector<double> &mz0_out)
{
  //This function reads initialises the macrocell spin configuration from a file.
  std::string cfgFile="initial_configuration.txt";
  std::ifstream f1(cfgFile);

  if(!f1)
  {
        std::cout << std::endl << "Failed to open the file: "<< cfgFile <<std::endl;
        exit(1);
  }

  for(int cell=0; cell<n_cells; cell++)
  {
    f1 >> mx0_out[cell];
    f1 >> my0_out[cell];
    f1 >> mz0_out[cell];

    std::cout<<"cell: "<<cell<<" mx0: "<<mx0_out[cell]<<" my0: "<<my0_out[cell]<<" mz0: "<<mz0_out[cell]<<"\n";


  }



  return 0;
}

int antiparallel_config(int n_cells,
                       std::vector<double> mx0_in,std::vector<double> my0_in,std::vector<double> mz0_in,
                       std::vector<double> &mx0_out,std::vector<double> &my0_out,std::vector<double> &mz0_out,
                       std::vector<int>material_id)
{
  //This function will be called when I need to create a DW. It will initialise the macrospin system in
  //an antiparallel config. w.r.t the middle region.

  int mat;
  for(int cell=0; cell<n_cells/2; cell++)
  {
    mat=material_id[cell];
    mx0_out[cell]=mx0_in[mat];
    my0_out[cell]=my0_in[mat];
    mz0_out[cell]=mz0_in[mat];

    std::cout<<"cell: "<<cell<<" mx0: "<<mx0_out[cell]<<" my0: "<<my0_out[cell]<<" mz0: "<<mz0_out[cell]<<"\n";
  }

  for(int cell=n_cells/2; cell<n_cells; cell++)
  {
    mat=material_id[cell];
    mx0_out[cell]=mx0_in[mat];
    my0_out[cell]=my0_in[mat];
    mz0_out[cell]=-mz0_in[mat];

    std::cout<<"cell: "<<cell<<" mx0: "<<mx0_out[cell]<<" my0: "<<my0_out[cell]<<" mz0: "<<mz0_out[cell]<<"\n";
  }

  return 0;
}

int uniform_config(int n_cells,
                       std::vector<double> mx0_in,std::vector<double> my0_in,std::vector<double> mz0_in,
                       std::vector<double> &mx0_out,std::vector<double> &my0_out,std::vector<double> &mz0_out,
                       std::vector<int>material_id)
{ //This function sets the initial configuration of the spin vectors mx0_out my0_out mz0_out from the macrospin file
  //using the mx_0 my_0 mz_0 vectors  from the input file

  int mat;
      
  for(int cell=0; cell<n_cells; cell++)
  {
    mat=material_id[cell];
    mx0_out[cell]=mx0_in[mat];
    my0_out[cell]=my0_in[mat];
    mz0_out[cell]=mz0_in[mat];

    //std::cout<<"cell: "<<cell<<" mx0: "<<mx0_out[cell]<<" my0: "<<my0_out[cell]<<" mz0: "<<mz0_out[cell]<<"\n";
  }
  
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


std::vector<int> interaction_list;
std::vector<int> start_neighbours;
std::vector<int> end_neighbours;
std::vector<int> count_int;
std::vector<int> left_edge_spins;
std::vector<int> right_edge_spins;

//---Functions
int check_edge_spins(int n_cells, 
                     std::vector<unsigned long long int>nx_cells_mat,
                     std::vector<int>material_id,
                     std::vector<int> count_int,
                     std::vector<double>xcoord,
                     std::vector<double>ycoord,
                     std::vector<double>zcoord,
                     std::vector<int> &left_edge_spins,
                     std::vector<int> &right_edge_spins)
{
  //This function checks which are the left and right edge spins and saves their ids in two int vectors
  double rtol = 1e-5;
  int mat_id;
  int nx_cells;


  for(int cell=0; cell<n_cells; cell++)
  {
    mat_id=material_id[cell]; //Get material id of cell
    nx_cells=nx_cells_mat[mat_id]; //Get number of nx cells for this material

    if(count_int[cell] <=2) //1 for 1 d, 2 for 2d , 3 for 3d 
    {

      if(fabs(xcoord[cell]-0.0)<rtol) //if the spin's coordinate corresponds the left edge
      {
        left_edge_spins.push_back(cell);
      }

      if(fabs(xcoord[cell]-(nx_cells-1.0))<rtol) //if the spin coordinate corresponds to the right edge..
      {

        right_edge_spins.push_back(cell);

      }

    }


  }



  return 0;
}


int create_interaction_list(int n_cells,
                            std::vector<double> xcoord,
                            std::vector<double> ycoord,
                            std::vector<double> zcoord,
                            std::vector<int> &int_list,
                            std::vector<int> &start,
                            std::vector<int> &end,
                            std::vector<int> &count_int)
{ //This function creates the interaction list.
  //The start and end vectors will help me find the neighbours of each macrospin cell in the interaction list.
  start.resize(n_cells);
  end.resize(n_cells);
  count_int.resize(n_cells);

  double r0=1.1;//range of interactions. This is currently equal to 1.1,meaning only nearest neighbours are counted. 
                //Atomic positions are incremented using a step of 1 in vectors xcoord, ycood, zcoord
  const double rtoll=r0*1.0e-5; //range tolerance

  for(int i=0; i<n_cells; i++) //For each macrospin cell
  {
   int count_interactions=0; //Set an interactions counter for each new macrospin cell
   for(int j=0; j<n_cells; j++) //Loop all the other cells
    {
      
      if(i==j)continue;//exclude self interactions
      else
      {
      //compute the distance between macrocells
      double rij=sqrt( pow((xcoord[i]-xcoord[j]),2.0) +
                       pow((ycoord[i]-ycoord[j]),2.0) +
                       pow((zcoord[i]-zcoord[j]),2.0) );
      if(rij-r0<rtoll)//if the distance agrees with the interaction range..
      {
       count_interactions++; //Count it 
       int_list.push_back(j);//Update interaction list with neighbour id
       int index_j=std::distance(int_list.begin(),int_list.end());
       end[i]=index_j;//Update end list
      }
      }
    }
    start[i]=end[i]-count_interactions; //Update start list
    count_int[i]=count_interactions; //count the number of interactions for each cell
  }


  return 0;
}

int generate_crystal_structure_f(int n_materials,
                                 int n_cells,
                                 std::vector<unsigned long long int>nx, std::vector<unsigned long long int>ny, std::vector<unsigned long long int>nz,
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
       //std::cout<<"Cell:"<<count_cell<<" Xcoord:"<<x[count_cell]<<" Ycoord:"<<y[count_cell]<<" Zcoord:"<<z[count_cell]<<" Material ID:"<<mat_id[count_cell]<<"\n";

       
       count_cell++;
      
      }
     }
    }
    base += nz[mat];
  }
  return 0;  
}

}

namespace macrospin{


  namespace internal{

    int alloc_memory(int n_cells)
    {


      macrospin::mx_0.resize(n_cells);
      macrospin::my_0.resize(n_cells);
      macrospin::mz_0.resize(n_cells);

      macrospin::mx.resize(n_cells);
      macrospin::my.resize(n_cells);
      macrospin::mz.resize(n_cells);

      return 0;
    }

    int set_initial_config(int n_cells)
    {
      if(input::force_DW_formation==0 && input::read_config_from_file==0)
      {
        macrospin::uniform_config(n_cells,
                                  input::mx_0, input::my_0, input::mz_0,
                                  macrospin::mx_0, macrospin::my_0, macrospin::mz_0,
                                  material::id);
      }
      else
      {  
        
        if(input::force_DW_formation==1)
        {
          macrospin::antiparallel_config(n_cells,
                                         input::mx_0, input::my_0, input::mz_0,
                                         macrospin::mx_0, macrospin::my_0, macrospin::mz_0,
                                         material::id);
        }
        
        if(input::read_config_from_file==1)
        {
          macrospin::read_config_from_file(input::n_cells,
                                           macrospin::mx_0, macrospin::my_0, macrospin::mz_0);
        }
      }
      return 0;
    }
  }
}

//---End of structure.cpp file.

