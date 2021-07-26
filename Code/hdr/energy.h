//---This is the energy.h file.
#pragma once

#ifndef ENERGY_H
#define ENERGY_H

//---Standard libraries


//---User-defined libraries


//---Namespace energy

namespace energy{
	
//---Variables
extern double uniax_anis;
extern double zeeman;

//---Functions
double uniax_anis_f(double mx, double my, double mz,
				    double ex, double ey, double ez,
				    double K,double V);


double zeeman_f(double mx, double my, double mz,
			    double bx, double by, double bz,
			    double B, double Ms, double V);

}

#endif
//---End of energy.h file.