//---This is the energy.cpp file. It defines the energy.h functions.


//---Standard libraries
#include<cmath>

//---User-defined libraries
#include"energy.h"
#include"input.h"

//---Namespace energy

namespace energy{
	
//---Variables
double uniax_anis=0.0;
double zeeman=0.0;


//---Functions 
 double uniax_anis_f(double mx, double my, double mz,
 				     double ex, double ey, double ez,
				     double K, double V)
 {
	return -K*V*pow(mx*ex + my*ey + mz*ez,2.0); //Uni-axial magnetocrystalline anisotropy energy: [J]
 }


 double zeeman_f(double mx, double my, double mz,
 			     double bx, double by, double bz,
			     double B, double Ms, double V)
 {


 	return -V*Ms*B*(mx*bx + my*by + mz*bz); //Zeeman energy: [J]

 }


}


//---End of energy.cpp file.

