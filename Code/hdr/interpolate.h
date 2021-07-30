//---This is the interpolate.h file.
#pragma once

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

//---Standard libraries
#include<vector>

//---User-defined libraries


//---Namespace cubicspline

namespace cubicspline{
	
//---Variables
extern std::vector<std::vector<double>>x_interpol,y_interpol;
extern std::vector<std::vector<double>>b,c,d;

extern std::vector<double>x_temp, y_temp;
extern std::vector<double>b_temp,c_temp,d_temp;

//---Functions
double polynome(double a, double b, double c,double d, double x, double xj);
int get_xinterval_id(std::vector<double>x_arr, double x, int &id);
int interpolate(int n,std::vector<double> x,std::vector<double> a,
                std::vector<double> &b, std::vector<double> &c, std::vector<double> &d);


}

#endif
//---End of interpolate.h file.