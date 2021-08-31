//---This is the utils.cpp file. It defines the utils.h functions.


//---Standard libraries
#include<vector>
#include<cmath>

//---User-defined libraries
#include"utils.h"

//---Namespace utils

namespace utils{
	
//---Variables



//---Functions

int max_element_1D_vec(std::vector<double> vec, double &maxx, int &id, bool absolute)
{	//finds the maximum element of vector // also returns the id of the maximum element
    //the absolute bool controls whether the elements of the vector will be considered in absolute value or not
	
	if(absolute==true)vec[0]=fabs(vec[0]);
	maxx=vec[0];
	id=0;

	for(int i=1; i<vec.size(); i++)
	{	
		if(absolute==true)vec[i]=fabs(vec[i]);
		if(vec[i]>maxx){
						maxx=vec[i];
						id=i;
					   }
	}

	return 0;
}

}


//---End of utils.cpp file.

