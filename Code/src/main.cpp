//---This is the main.cpp file. It calls the main routines of the code.


//---Standard libraries
#include<iostream>

//---User-defined libraries
#include"input.h"
#include"solver.h"
#include"particle.h"
#include"field.h"
#include"init.h"
#include"tempscaling.h"
#include"interpolate.h"
#include"output.h"
#include"simulation.h"

int main()
{


init::internal::run();

init::internal::end();

return 0;
}





//---End of main.cpp file
