//---This is the rootfind.h file.
#pragma once

#ifndef ROOTFIND_H
#define ROOTFIND_H

//---Standard libraries


//---User-defined libraries


//---Namespace rootfind

namespace rootfind{
	
//---Variables
extern double x0;
extern double N_iter;
extern double TOL;

//---Functions
double NewtonRaphson(double (*f)(double, double, double, double),
	                 double (*g)(double, double, double, double),
	                 double x0, double TOL, int N_iter,
	                 double eps, 
	                 double T, double Tc);


}

#endif
//---End of rootfind.h file.