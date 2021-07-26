//---This is the field.h file.
#pragma once

#ifndef FIELD_H
#define FIELD_H

//---Standard libraries


//---User-defined libraries


//---Namespace field

namespace field{
	
//---Variables
extern double Bx_ani, By_ani, Bz_ani;
extern double Bx_app, By_app, Bz_app;
extern double Bx_lon, By_lon, Bz_lon;
extern double Bx_eff, By_eff, Bz_eff;

extern double torque_x, torque_y, torque_z, torque_mod;

//---Functions
int uniax_anis_f(double mx, double my, double mz,
			     double ex, double ey, double ez,
				 double K, double Ms,
				 double &Bx_ani, double &By_ani, double &Bz_ani);
				

int zeeman_f(double B_app, double bx, double by, double bz,
		     double &Bx_app, double &By_app, double &Bz_app);


int longitudinal_f(double mx, double my, double mz,
				   double chi_par, double m_e,
				   double &Bx_lon, double &By_lon, double &Bz_lon);

int effective_f(double Bx_ani, double By_ani, double Bz_ani,
				double Bx_app, double By_app, double Bz_app,
				double &Bx_eff, double &By_eff, double &Bz_eff);

int effective_torque(double mx, double my, double mz,
 					 double Bx_eff, double By_eff, double Bz_eff,
					 double &torque_x, double &torque_y, double &torque_z, double  &torque_mod);

int calculate();


}

#endif
//---End of field.h file.