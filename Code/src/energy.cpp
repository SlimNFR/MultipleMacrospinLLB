//---This is the energy.cpp file. It defines the energy.h functions.


//---Standard libraries
#include<cmath>
#include<vector>

//---User-defined libraries
#include"energy.h"
#include"input.h"
#include"structure.h"

//---Namespace energy

namespace energy{
	
//---Variables

//Total energies
double anis_total=0.0;
double exchange_total=0.0;
double zeeman_total=0.0;
double total=0.0;


//Energies per cell
std::vector<double> anis_cell;
std::vector<double> exchange_cell;
std::vector<double> zeeman_cell;
std::vector<double> total_cell;

//---Functions 
 int anis_cell_f(int n_cells,
                 std::vector<double> K, std::vector<double> V,
                 std::vector<double> ex, std::vector<double> ey, std::vector<double> ez,
                 std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                 std::vector<int> mat_id,
                 std::vector<double> &anis_cell)
 {
    int mat; //stores the material id of the cell

    for(int cell=0; cell<n_cells; cell++)
    {//loop cells
        mat=mat_id[cell];
        anis_cell[cell] = -K[mat]*V[mat]*pow((mx[cell]*ex[mat] + my[cell]*ey[mat] + mz[cell]*ez[mat]),2.0); //E_ani : [J]

    }
    
    return 0; 
 }


 int zeeman_cell_f(int n_cells,
                   double B, double bx, double by, double bz,
                   std::vector<double> Ms, std::vector<double> V,
                   std::vector<double> mx, std::vector<double> my, std::vector<double> mz,
                   std::vector<int> mat_id,
                   std::vector<double> &zeeman_cell)
 {
   int mat; //stores the material id of the cell

   for(int cell=0; cell<n_cells; cell++)
   {

      mat=mat_id[cell];
      zeeman_cell[cell] = -V[mat]*Ms[mat]*B*(mx[cell]*bx + my[cell]*by + mz[cell]*bz); //E_zeeman : [J]

   }

 	return 0;

 }


 int exchange_cell_f(int n_cells,double lengthscale,
                     std::vector<unsigned long long int> macrocell_size,std::vector<std::vector<double>>A_T_matrix,
                     std::vector<int> int_list, std::vector<int> start_neighbours, std::vector<int> end_neighbours,
                     std::vector<int> material_id,
                     std::vector<double> mx, std::vector<double> my,std::vector<double> mz,
                     std::vector<double>& exc_cell)
 {
   //beware, currently, the exchange energy is implemented for a perfect cuboid, hx=hy=hz.
   
   std::fill(exc_cell.begin(), exc_cell.end(), 0.0);

   //Temp. variables
   int neighbour, mat_id_neighbour, mat_id_cell;
   double A;


   //Loop each cell
   for(int cell=0; cell<n_cells; cell++)
   {

      mat_id_cell=material_id[cell]; //Get material ID of cell
      
      //Count all cell's neighbours
      for(int count_neighbour=start_neighbours[cell]; count_neighbour<end_neighbours[cell]; count_neighbour++) 
      {

         neighbour = int_list[count_neighbour]; //Get neighbour from interaction list
         mat_id_neighbour=material_id[neighbour]; //Get material id of neighbour

         A = A_T_matrix[mat_id_cell][mat_id_neighbour]; //Get exchange from exchange matrix
         
         //Calculate exchange energy

         double dprod_cell_neigh=(mx[neighbour]*mx[cell] + my[neighbour]*my[cell] + mz[neighbour]*mz[cell]);
         double dprod_cell_cell=(mx[cell]*mx[cell] + my[cell]*my[cell] + mz[cell]*mz[cell]);

         exc_cell[cell] = - 2.0*A*macrocell_size[mat_id_cell]*lengthscale*(dprod_cell_neigh - dprod_cell_cell); // [Joules]
         
      }

   }

   return 0;
 }

 int anis_total_f(int n_cells,
                  std::vector<double> anis_cell,
                  double &anis_total)
 {
   //calculates the total anisotropy energy

   anis_total = 0.0;
   for(int cell=0; cell<n_cells; cell++)
   {
      anis_total +=anis_cell[cell];
   }

   return 0;
 }

int zeeman_total_f(int n_cells,
                   std::vector<double> zeeman_cell,
                   double &zeeman_total)
 {
   //calculates the total zeeman energy

   zeeman_total = 0.0;
   for(int cell=0; cell<n_cells; cell++)
   {
      zeeman_total +=zeeman_cell[cell];
   }

   return 0;
 }


 int exchange_total_f(int n_cells,
                      std::vector<double> exc_cell,
                      double &exc_total)

 { //calculates the total exchange energy

   exc_total = 0.0;
   for(int cell=0; cell<n_cells; cell++)
   {
      exc_total +=exc_cell[cell];
   }

   return 0;
 }

 int total_f(double ani,
             double exc,
             double zee,
             double & total)
 {

   total = ani + exc + zee;

   return 1;

 }

 namespace internal {


      int alloc_memory(int n_cells)
      { //This function allocates memory for the

         energy::exchange_cell.resize(n_cells,0.0);
         energy::zeeman_cell.resize(n_cells,0.0);
         energy::anis_cell.resize(n_cells,0.0);
         

         return 0;
      }

      int calculate()
      {

         energy::anis_cell_f(input::n_cells,
                             input::K_T, input::macrocell_volume,
                             input::ex, input::ey, input::ez,
                             macrospin::mx, macrospin::my,macrospin::mz,
                             material::id,
                             energy::anis_cell);

         energy::zeeman_cell_f(input::n_cells,
                               input::B_app, input::bx, input::by, input::bz, 
                               input::Ms_T, input::macrocell_volume,
                               macrospin::mx, macrospin::my, macrospin::mz,
                               material::id,
                               energy::zeeman_cell);

         energy::exchange_cell_f(input::n_cells, input::lengthscale,
                                 input::macrocell_size,input::A_T_matrix, 
                                 material::interaction_list, material::start_neighbours, material::end_neighbours, material::id,                                 
                                 macrospin::mx, macrospin::my, macrospin::mz, 
                                 energy::exchange_cell);

         energy::anis_total_f(input::n_cells,
                              energy::anis_cell,
                              energy::anis_total);

         energy::zeeman_total_f(input::n_cells,
                                energy::zeeman_cell,
                                energy::zeeman_total);

         energy::exchange_total_f(input::n_cells,
                                  energy::exchange_cell,
                                  energy::exchange_total);

         energy::total_f(energy::anis_total,
                         energy::exchange_total,
                         energy::zeeman_total,
                         energy::total);




         return 0;
      }


 }


}


//---End of energy.cpp file.

