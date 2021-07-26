//---This is the rootfind.cpp file. It defines the rootfind.h functions.


//---Standard libraries
#include<cmath>

//---User-defined libraries
#include"rootfind.h"

//---Namespace rootfind

namespace rootfind{
	
//---Variables
double x0 = 1000.0; //initial guess value
double N_iter=200;  //allowed number of iterations for the root find algorithm
double TOL=1e-10;   //if f(x)<TOL then x is solution.

//---Functions
double NewtonRaphson(double (*f)(double, double, double, double),
                                 double (*g)(double, double, double, double),
                                 double x0, double TOL, int N_iter,
                                 double eps,
                                 double T, double Tc)
{//Root-finding algorithm based on the Newton-Raphson method
 //Algorithm: xn=x0-f(x0)/g(x0) where f(x0) is the function evaluated at x0 and g(x0) is the first derivative evaluated at x0

        double xn;

        if(f(x0,T,Tc,eps) != 0.0) //If the initial point is not solution, then search
    {
        for(int i=1;i<=N_iter;i++)
        {
            xn=x0-(f(x0,T,Tc,eps)/g(x0,T,Tc,eps));
            if(fabs(f(xn,T,Tc,eps))<TOL)return (xn*eps*T)/(3*Tc); //If xn is solution then stop 
            else
            {
                //d=fabs((xn-x0)/xn)*100 %; Error (% percentage)
                x0=xn; //the new initial point is swapped with b
            }
        }
        //if for loop didn't break, consider the solution to be the last point.
        return (xn*eps*T)/(3*Tc);

    }
    else{
        //
        return (x0*eps*T)/(3*Tc);
    }

return 1;//if it didn't work return 1

}



}


//---End of rootfind.cpp file.

