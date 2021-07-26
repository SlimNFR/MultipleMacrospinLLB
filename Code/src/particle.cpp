//---This is the particle.cpp file. It defines the particle.h functions.


//---Standard libraries
#include<vector>

//---User-defined libraries
#include"particle.h"
//---Namespace particle

namespace particle{
	
//---Variables
double mx;
double my;
double mz;

double checkpoint_mx;
double checkpoint_my;
double checkpoint_mz;


//---Functions
int update_magnetisation(double mx_in,double my_in, double mz_in,
						 double &mx, double &my, double &mz)
{
	mx = mx_in; my = my_in; mz = mz_in;
	return 0;
}

}


//---End of particle.cpp file.

