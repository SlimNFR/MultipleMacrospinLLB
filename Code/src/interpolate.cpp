//---This is the interpolate.cpp file. It defines the interpolate.h functions.


//---Standard libraries
#include<vector>
#include<cmath>

//---User-defined libraries
#include"interpolate.h"

//---Namespace cubicspline

namespace cubicspline{
	
//---Variables
std::vector<double>x_interpol,y_interpol;
std::vector<double>b,c,d;


//---Functions
double polynome(double a, double b, double c,double d, double x, double xj)
{	//This is the 3rd degree polynome used in the Cubic Spline Interpolation method
    return a + b*(x-xj) +c*pow((x-xj),2.0) + d*pow((x-xj),3.0);
}


int get_xinterval_id(std::vector<double>x_arr, double x, int &id)
{	//This function determines the spline interval corresponding to the x variable. 

    for(long unsigned int i=0;i<x_arr.size()-1;i++)
    {

        if(x_arr[i]<= x && x <= x_arr[i+1])
        {

            id = i;
            return 0;
        }


    }
    return 0;

}


int interpolate(int n, std::vector<double> x,std::vector<double> a,
                std::vector<double> &b, std::vector<double> &c, std::vector<double> &d)
{	//This function interpolates the given x and a set of points using the Cubic Spline method
    //Array named a is y which in turn represents f evaluated at x!
    /** Numerical Analysis 9th ed - Burden, Faires (Ch. 3 Natural Cubic Spline, Pg. 149) */

    std::vector<double> h(n), A(n), l(n + 1),u(n + 1), z(n+1);

    // Step 1 /
    for (int i = 0; i <= n - 1; ++i) h[i] = x[i + 1] - x[i];

    // Step 2 /

    for (int i = 1; i <= n - 1; ++i)
        A[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];

    // Step 3 //
    l[0] = 1;
    u[0] = 0;
    z[0] = 0;

    // Step 4 /
    for (int i = 1; i <= n - 1; ++i) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    // Step 5 /
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    // Step 6 /
    for (int j = n - 1; j >= 0; --j) {
        c[j] = z[j] - u[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    return 0;
}


}


//---End of interpolate.cpp file.

