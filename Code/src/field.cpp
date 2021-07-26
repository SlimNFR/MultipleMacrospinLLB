//---This is the field.cpp file. It defines the field.h functions.


//---Standard libraries
#include<iostream>
#include<cmath>

//---User-defined libraries
#include"field.h"
#include"particle.h"
#include"input.h"

//---Namespace field

namespace field{
	
//---Variables
double Bx_ani, By_ani, Bz_ani = 0.0;
double Bx_app, By_app, Bz_app = 0.0;
double Bx_lon, By_lon, Bz_lon = 0.0;
double Bx_eff, By_eff, Bz_eff = 0.0;

double torque_x, torque_y, torque_z, torque_mod = 0.0;

//---Functions
int uniax_anis_f(double mx, double my, double mz,
				 double ex, double ey, double ez,
				 double K, double Ms,
				 double &Bx_ani, double &By_ani, double &Bz_ani)
{
	//calculates the anisotropy field components
	
	double Hk=2.0*K/Ms; //Anisotropy field strength [T]

	//std::cout<<"Ms inside uniax field func: "<<Ms<<"\n";
	//std::cout<<"Withing anis field func: | mx:"<<mx<<"| my:"<<my<<"| mz:"<<mz<<"\n";
	Bx_ani = Hk*(mx*ex + my*ey + mz*ez)*ex;
	By_ani = Hk*(mx*ex + my*ey + mz*ez)*ey;
	Bz_ani = Hk*(mx*ex + my*ey + mz*ez)*ez;

	//std::cout<<"Within function"<<Bx_ani<<" "<<By_ani<<" "<<Bz_ani<<"\n";
	//std::cout<<"Within function"<<mx<<" "<<my<<" "<<mz<<"\n";
	//std::cout<<"Within function"<<ex<<" "<<ey<<" "<<ez<<"\n";

	return 0;
}
				

int zeeman_f(double B_app, double bx, double by, double bz,
			 double &Bx_app, double &By_app, double &Bz_app)
{
	//calculates the zeeman field components
	//std::cout<<"Inside zeeman field function"<<"bx: "<<bx<<"|by: "<<by<<"|bz: "<<bz<<"|B_app: "<<B_app<<"\n";
	Bx_app=B_app*bx;
	By_app=B_app*by;
	Bz_app=B_app*bz;

	return 0;
}


int longitudinal_f(double mx, double my, double mz,
				   double chi_par, double m_e,
				   double &Bx_lon, double &By_lon, double &Bz_lon)
{
	//calculates the longitudinal field components
	//I assume I work below Tc for now

	if(input::T==0){Bx_lon=0;By_lon=0;Bz_lon=0; return 0;}
	//Temporary variables
	double m_squared = mx*mx + my*my + mz*mz;
	//std::cout<<"Withing long field func: | mx:"<<mx<<"| my:"<<my<<"| mz:"<<mz<<"\n";
	//chi_par=0.1;
	double pre_factor = (0.5*(1.0/chi_par)) * (1.0 - (m_squared/(m_e*m_e)) );
	//std::cout<<"Withing longitudinal field:"<<"chi_par: "<<chi_par<<"| pre_factor: "<<pre_factor<<"| m_squared: "<<m_squared<<"| m_e:"<<m_e<<"\n";
	Bx_lon = pre_factor*mx;
	By_lon = pre_factor*my;
	Bz_lon = pre_factor*mz;

	return 0;

}


int effective_f(double Bx_ani, double By_ani, double Bz_ani,
				double Bx_app, double By_app, double Bz_app,
				double &Bx_eff, double &By_eff, double &Bz_eff)
{
	//calculates the total field (effective) components

	Bx_eff = Bx_ani + Bx_app + Bx_lon;
	By_eff = By_ani + By_app + By_lon;
	Bz_eff = Bz_ani + Bz_app + Bz_lon;

	return 0;

}

int effective_torque(double mx, double my, double mz,
					 double Bx_eff, double By_eff, double Bz_eff,
					 double &torque_x, double &torque_y, double &torque_z, double  &torque_mod)
{

	torque_x = (my*Bz_eff - mz*By_eff);
	torque_y = (mz*Bx_eff - mx*Bz_eff);
	torque_z = (mx*By_eff - my*Bx_eff);

	torque_mod = sqrt(torque_x*torque_x+
				 	  torque_y*torque_y+
				 	  torque_z*torque_z);

	return 0;
}

int calculate()
{
	field::uniax_anis_f(particle::mx, particle::my, particle::mz,
						input::ex, input::ey, input::ez,
						input::K_T, input::Ms_T, 
						field::Bx_ani, field::By_ani, field::Bz_ani);

	field::zeeman_f(input::B_app, input::bx, input::by, input::bz,
					field::Bx_app, field::By_app, field::Bz_app);

	field::longitudinal_f(particle::mx, particle::my, particle::mz,
			   		      input::chi_par, input::m_e,
				          field::Bx_lon, field::By_lon, field::Bz_lon);

	field::effective_f(field::Bx_ani, field::By_ani, field::Bz_ani,
					   field::Bx_app, field::By_app, field::Bz_app,
					   field::Bx_eff, field::By_eff, field::Bz_eff);

	field::effective_torque(particle::mx, particle::my, particle::mz,
							field::Bx_eff, field::By_eff, field::Bz_eff,
					 		field::torque_x, field::torque_y, field::torque_z, field::torque_mod);

	//std::cout<<"Bz_ani:"<<Bz_ani<<"| Bz_app"<<" "<<Bz_app<<"|Bz_lon: "<<Bz_lon<<"\n";

	return 0;
}

}


//---End of field.cpp file.

