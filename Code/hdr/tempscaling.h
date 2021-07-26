//---This is the tempscaling.h file.
#pragma once

#ifndef TEMPSCALING_H
#define TEMPSCALING_H

//---Standard libraries
#include<iostream>
#include<fstream>
#include<vector>

//---User-defined libraries


//---Namespace tempscaling

namespace tempscaling{
	
//---Variables


//---Functions
//obtain instantaneous parameters
int alpha_par_f(double T, double Tc, double lambda,double &alpha_par);
int alpha_perp_f(double T, double Tc, double lambda, double &alpha_perp);
int equilibrium_magn_f(double T, std::vector<double> x_interpol, std::vector<double> y_interpol,
					   std::vector<double> b,std::vector<double> c,std::vector<double> d, double &m_e);
int chi_par_f(double (*dLangevin)(double, double, double, double),
			  double T, double Tc, double mu_s, double m_e, double eps, double &chi_par);

int K_at_T_f(double K0_SI, double m_e, double &K_T);
int Ms_at_T_f(double Ms0_SI, double m_e, double &Ms_T);


//obtain equilibrium curves
int K_vs_T_curve_f(double Tc, double K0_SI,
				   std::vector<double> x_interpol, std::vector<double>y_interpol,
				   std::vector<double> b, std::vector<double> c, std::vector<double> d,std::ofstream &f1);

int chipar_vs_T_curve_f(double Tc, double eps, 
						double mu_s,
						std::vector<double>x_interpol, std::vector<double>y_interpol,
						std::vector<double>b, std::vector<double>c, std::vector<double>d,
						std::ofstream &f1);

int m_vs_T_curve_f(double Tc,
				   std::vector<double>x_interpol, std::vector<double>y_interpol,
				   std::vector<double>b, std::vector<double>c, std::vector<double>d,			 	
				   std::ofstream &f1, std::ofstream &f2);

//others..
int get_mVsT_points_to_interpol(double Tc, double eps, std::vector<double> &x_interpol, std::vector<double> &y_interpol);
int interpolate_mVsT_points(std::vector<double>x_interpol, std::vector<double>y_interpol,
							std::vector<double> &b, std::vector<double> &c, std::vector<double> &d);

}


namespace tempscaling{



	namespace internal{
		int obtain_interpolation_polynome_mVsT_data();
		int call_mVsT_sim();
		int call_chiparVsT_sim();
		int call_KVsT_sim();
		int calc_parameters_at_T();
	}
}

#endif
//---End of tempscaling.h file.
