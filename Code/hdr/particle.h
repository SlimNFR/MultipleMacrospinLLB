//---This is the particle.h file.
#pragma once

#ifndef PARTICLE_H
#define PARTICLE_H

//---Standard libraries
#include<vector>

//---User-defined libraries


//---Namespace particle

namespace particle{
    
//---Variables
extern double mx;
extern double my;
extern double mz;

extern double checkpoint_mx;
extern double checkpoint_my;
extern double checkpoint_mz;

//---Functions
int update_magnetisation(double mx_in,double my_in, double mz_in,
						 double &mx, double &my, double &mz);

}

#endif
//---End of particle.h file.