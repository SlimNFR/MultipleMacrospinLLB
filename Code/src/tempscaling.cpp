//---This is the tempscaling.cpp file. It defines the tempscaling.h functions.


//---Standard libraries
#include<fstream>
#include<iostream>

//---User-defined libraries
#include"tempscaling.h"
#include"rootfind.h"
#include"interpolate.h"
#include"equation.h"
#include"input.h"
#include"output.h"

//---Namespace tempscaling

namespace tempscaling{
	
//---Variables


//---Functions

int alpha_par_f(double T, double Tc, double lambda,double &alpha_par)
{	//calculates alpha_par at given temperature T

	alpha_par = (2.0/3.0)*(T/Tc)*lambda;
	return 0;

}
int alpha_perp_f(double T, double Tc, double lambda, double &alpha_perp)
{	//calculates alpha_perp at given temperature T
	
	alpha_perp = lambda*(1-(T/(3*Tc)));
	return 0;
}


int equilibrium_magn_f(double T, std::vector<double> x_interpol, std::vector<double> y_interpol,
					   std::vector<double> b,std::vector<double> c,std::vector<double> d, double &m_e)
{	//calculates the equilibrium magnetisation for a given temperature T

	int id;

	cubicspline::get_xinterval_id(x_interpol, T, id);
	m_e=cubicspline::polynome(y_interpol[id],b[id],c[id],d[id],T,x_interpol[id]);
	return 0;
}




int chi_par_f(double (*dLangevin)(double, double, double, double),
			  double T, double Tc, double mu_s, double m_e, double eps, double &chi_par)
{	//calculates the parallel susceptibility for a given temperature T

	double kb=input::k_B;
	double dL=dLangevin(T, Tc, eps,m_e);
	if(T==0)return 0;
	chi_par = (mu_s/kb)*(dL/(T-3*(Tc/eps)*dL));

	return 0;
}

int K_at_T_f(double K0_SI, double m_e, double &K_T)
{

	K_T = K0_SI * pow(m_e,3.0);

	return 0;
}


int Ms_at_T_f(double Ms0_SI, double m_e, double &Ms_T)
{

	Ms_T = Ms0_SI*m_e;

	return 0;
}


int K_vs_T_curve_f(double Tc, double K0_SI,
				   std::vector<double> x_interpol, std::vector<double>y_interpol,
				   std::vector<double> b, std::vector<double> c, std::vector<double> d,std::ofstream &f1)
{
	int T;
	double m_e;
	double K_T;

	for(T=0; T<=Tc; T++)
	{
		tempscaling::equilibrium_magn_f((double)T, x_interpol,y_interpol,b,c,d,m_e);
		tempscaling::K_at_T_f(K0_SI, m_e, K_T);
		f1<<T<<" "<<K_T/K0_SI<<"\n";

	}
	return 0;
}


int chipar_vs_T_curve_f(double Tc, double eps, 
						double mu_s,
						std::vector<double>x_interpol, std::vector<double>y_interpol,
						std::vector<double>b, std::vector<double>c, std::vector<double>d,
						std::ofstream &f1)
{
	int T;
	double m_e;
	double chi_par;

	//Loop temperatures and print to file
	for(T=0; T<=Tc; T++)
	{
		tempscaling::equilibrium_magn_f((double)T, x_interpol,y_interpol,b,c,d,m_e);
		tempscaling::chi_par_f(equation::Langevin_df, T, Tc, mu_s,m_e, eps, chi_par);
		f1<<T<<" "<<chi_par<<"\n";
	}

	return 0;
}


int m_vs_T_curve_f(double Tc, 
				   std::vector<double>x_interpol, std::vector<double>y_interpol,
				   std::vector<double>b, std::vector<double>c, std::vector<double>d,			 	
				   std::ofstream &f1, std::ofstream &f2)
{	//Calculates the whole m_e(T) curve.
	int T;

	//print to file the NR data
	for(long unsigned int i=0; i<x_interpol.size(); i++)
	{
		f1<<x_interpol[i]<<" "<<y_interpol[i]<<"\n";
	}

	int id; //This will store the interval id of each new x value (temperature value).
	for (T = 0; T <=Tc; T++)
	{	
	    cubicspline::get_xinterval_id(x_interpol, T, id);
	    //std::cout<<T<<" "<<cubicspline::y_interpol[id]<<" "<<cubicspline::b[id]<<" "<<cubicspline::c[id]<<" "<<cubicspline::d[id]<<"\n";

	    //print to file the interpolated data
	    f2<<T<<" "<<cubicspline::polynome(y_interpol[id],b[id],c[id],d[id],T,x_interpol[id])<<"\n"; //cubicspline::polynome(..)=m_e...
	}

	return 0;
}





int get_mVsT_points_to_interpol(double Tc, double eps, std::vector<double> &x_interpol, std::vector<double> &y_interpol)
{	//This function will obtain several m_e(T) points solving the Curie Weiss equation with the NR algorithm. 
	//These will later be interpolated.

	int T;
	double m_e;

	// Up to T=Tc-nearTc, I will solve the Curie Weiss equation in steps of broad_step_T.
	// For T>Tc - nearTc I use discrete_step_T.
	int broad_step_T=50;
	int discrete_step_T=1;
	int nearTc=50;


	//For my given temperature range, call the NR algorithm and solve the Curie-Weiss equation.
	for(T=0;T<=Tc;T++)
	{
	        if(T==0 || T==Tc)
	        {//The NR algorithm will not work for T=0 or T=Tc.

	                x_interpol.push_back(T);
	                if(T==0){m_e=1.0; y_interpol.push_back(m_e);}
	                if(T==Tc){m_e=0.0; y_interpol.push_back(m_e);}
	        }

	    else
	    {
	            if(T<=Tc-nearTc && T%broad_step_T == 0){

	                m_e=rootfind::NewtonRaphson(equation::CurieWeiss_f,equation::CurieWeiss_df,
	                							rootfind::x0,rootfind::TOL,rootfind::N_iter,eps,T,Tc);
	                x_interpol.push_back(T);
	                y_interpol.push_back(m_e);
	                //std::cout<<"x and y: "<<cubicspline::x_interpol.back()<<" "<<cubicspline::y_interpol.back()<<"\n";

	            }

	            else if (T>Tc-nearTc && T%discrete_step_T == 0){
	                m_e=rootfind::NewtonRaphson(equation::CurieWeiss_f,equation::CurieWeiss_df,
	                							rootfind::x0,rootfind::TOL,rootfind::N_iter,eps,T,Tc);
	                x_interpol.push_back(T);
	                y_interpol.push_back(m_e);
	                //std::cout<<"x and y: "<<cubicspline::x_interpol.back()<<" "<<cubicspline::y_interpol.back()<<"\n";
	            }
	        }

	}
	return 0;
}


int interpolate_mVsT_points(std::vector<double>x_interpol, std::vector<double>y_interpol,
							std::vector<double> &b, std::vector<double> &c, std::vector<double> &d)
{

//Interpolate the x_interpol,y_interpol data points via the Cubic spline algorithm to find the whole m_e(T) curve.

	int n = x_interpol.size(); //This will set the dimensions of my coefficient arrays b,c,d. 
	n--;
	b.resize(n);
	c.resize(n+1);
	d.resize(n);

	cubicspline::interpolate(n,x_interpol,y_interpol,b,c,d);

	return 0;
}


}

namespace tempscaling{



	namespace internal{


		int obtain_interpolation_polynome_mVsT_data()
		{   //this function will automatically set up the interpolated m_e(T) data
			tempscaling::get_mVsT_points_to_interpol(input::Tc, input::eps,
										 			 cubicspline::x_interpol, cubicspline::y_interpol);
			tempscaling::interpolate_mVsT_points(cubicspline::x_interpol, cubicspline::y_interpol,
												 cubicspline::b, cubicspline::c, cubicspline::d);

			return 0;
		}

		int call_mVsT_sim()
		{	//this function will run a simulation to obtain the full me_vs(T) curve.
			tempscaling::m_vs_T_curve_f(input::Tc,  
										cubicspline::x_interpol, cubicspline::y_interpol,
										cubicspline::b, cubicspline::c, cubicspline::d,
										output::file_Meq_temp_NR, output::file_Meq_temp_CS);
			return 0;
		}

		int call_chiparVsT_sim()
		{	//this function will run a simulation to obtain the full me_vs(T) curve.
			tempscaling::chipar_vs_T_curve_f(input::Tc, input::eps,input::mu_s,
											 cubicspline::x_interpol, cubicspline::y_interpol,
											 cubicspline::b,cubicspline::c, cubicspline::d, output::file_X_temp);
			return 0;
		}

		int call_KVsT_sim()
		{


			tempscaling::K_vs_T_curve_f(input::Tc, input::K0_SI, 
										cubicspline::x_interpol, cubicspline::y_interpol,
										cubicspline::b, cubicspline::c, cubicspline::d, output::file_K_temp);

			return 0;

		}


		int calc_parameters_at_T()
		{
			tempscaling::alpha_par_f(input::T, input::Tc, input::lambda, input::alpha_par);
			tempscaling::alpha_perp_f(input::T, input::Tc, input::lambda, input::alpha_perp);
			tempscaling::equilibrium_magn_f(input::T, cubicspline::x_interpol, cubicspline::y_interpol,
											cubicspline::b, cubicspline::c, cubicspline::d, input::m_e);
			tempscaling::chi_par_f(equation::Langevin_df, input::T, input::Tc, input::mu_s, input::m_e, input::eps, input::chi_par);
			tempscaling::Ms_at_T_f(input::Ms0_SI, input::m_e, input::Ms_T);
			tempscaling::K_at_T_f(input::K0_SI, input::m_e, input::K_T);
			
			return 0;

		}

	}

}



//---End of tempscaling.cpp file.

