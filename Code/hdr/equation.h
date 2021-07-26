//---This is the equation.h file.
#pragma once

#ifndef EQUATION_H
#define EQUATION_H

//---Standard libraries


//---User-defined libraries


//---Namespace equation

namespace equation{
	
//---Variables
extern double fx, fy, fz;//These variables hold the left-hand side value of the LLB equation on all 3 dimensions.

//---Functions
double CurieWeiss_f(double x, double T, double Tc, double eps);
double CurieWeiss_df(double x, double T, double Tc, double eps);
double Langevin_df(double T, double Tc, double eps, double m_e);
int LLB_classic(double mx, double my, double mz,
				double Bx_eff, double By_eff, double Bz_eff,
				double gamma, double alpha_par, double alpha_perp,
				double &fx,double &fy,double &fz);


}

#endif
//---End of equation.h file.