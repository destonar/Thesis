#ifndef MINKOPERATIONS_H
#define MINKOPERATIONS_H
#include <Eigen/Dense>
#include "ndim_const.h"
#include <set>
#include <vector>
#include <iostream>

using namespace Eigen;
using namespace std;

class MinkOperations
{
	Array<int, Dynamic, GLOBAL_NDIM> x1;
	Array<int, Dynamic, GLOBAL_NDIM> x2;
public:
	MinkOperations(const Array<int, Dynamic, GLOBAL_NDIM>& in_x1, const Array<int, Dynamic, GLOBAL_NDIM>& in_x2);

	Array<int, Dynamic, GLOBAL_NDIM> minksum();
};


#endif