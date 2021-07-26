//---This is the equation.cpp file. It defines the equation.h functions.


//---Standard libraries
#include<cmath>
#include<iostream>

//---User-defined libraries
#include"equation.h"



//---Namespace equation

namespace equation{
	
//---Variables

double fx, fy, fz = 0.0; //in these variables I will save the right hand side value of the LLB equation.

//---Functions
double CurieWeiss_f(double x, double T, double Tc, double eps)
{//This is the Curie Weiss self consistent equation
    return ( (eps/3.0)*(T/Tc)*x - 1.0/(tanh(x)) + 1.0/x );
}
double CurieWeiss_df(double x, double T, double Tc, double eps)
{//The derivative of the Curie Weiss equation
    return ( (eps/3.0)*(T/Tc)  + pow( (1.0/sinh(x)),2.0 ) - 1.0/(x*x)  );
}

double Langevin_df(double T, double Tc, double eps, double m_e)
{//The derivative of the Langevin function L(x)=coth(x)-1/x. L'(x)=-1/sinh(x)**2 +1/x**2
    double x ;
    if(T<1e-5)return 0; //if the temperature is close to 0 then the x argument shoots to +inf. limit of L'(x) will go to 0.
    if(T==Tc)return 1.0/3.0; //if T==Tc then L'(x) is 1/3
    else
    {
        x = (3.0/eps)*(Tc/T)*m_e;
        return ((-1.0/pow(sinh(x),2.0)) + 1.0/pow(x,2.0)) ;

    }
    

    
}

int LLB_classic(double mx, double my, double mz,
                double Bx_eff, double By_eff, double Bz_eff,
                double gamma, double alpha_par, double alpha_perp,
                double &fx,double &fy,double &fz)
{
    //sets up the LLB equation

    //Temporary variables
    double m_squared = mx*mx + my*my + mz*mz;
    double mTimesB = mx*Bx_eff + my*By_eff + mz*Bz_eff;
    double C1 = gamma*alpha_par/m_squared;
    double C2 = gamma*alpha_perp/m_squared;
    

    //std::cout<<"|gamma: "<<gamma<<"|alpha_par: "<<alpha_par<<"|alpha_perp: "<<alpha_perp<<"\n";
    // std::cout<<"Bx_eff: "<<Bx_eff<<"|By_eff: "<<By_eff<<"|Bz_eff: "<<Bz_eff
    //    <<"|mx: "<<mx<<"|my: "<<my<<"|mz: "<<mz<<"\n";

    fx = +gamma*(By_eff*mz - my*Bz_eff) + C1*(mx*Bx_eff)*mx + C2*(Bx_eff*m_squared - mx*mTimesB);
    fy = -gamma*(Bx_eff*mz - mx*Bz_eff) + C1*(my*By_eff)*my + C2*(By_eff*m_squared - my*mTimesB);
    fz = +gamma*(Bx_eff*my - mx*By_eff) + C1*(mz*Bz_eff)*mz + C2*(Bz_eff*m_squared - mz*mTimesB);

    
    //std::cout<<"Bx_eff: "<<Bx_eff<<"|By_eff: "<<By_eff<<"|Bz_eff: "<<Bz_eff
    //    <<"|mx: "<<mx<<"|my: "<<my<<"|mz: "<<mz<<"\n";
    
    //std::cout<<"|fx: "<<fx<<"|fy: "<<fy<<"|fz: "<<fz<<"\n";

    return 0;
}


}


//---End of equation.cpp file.

