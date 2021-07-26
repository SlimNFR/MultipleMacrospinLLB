//---This is the init.h file.
#pragma once

#ifndef INIT_H
#define INIT_H

//---Standard libraries
#include<chrono>

//---User-defined libraries
#include"tempscaling.h"

//---Namespace init

namespace init{
	
//---Variables
extern std::chrono::high_resolution_clock::time_point RUN_TIME_START;
extern std::chrono::high_resolution_clock::time_point RUN_TIME_END;
extern double RUN_TIME_TOTAL;

//---Functions
int parameters();
int files();
int sim();

}

namespace init{


	namespace internal{
		int run();
		int end();
	}
}

#endif
//---End of init.h file.
