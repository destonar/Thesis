#ifndef GRID_H
#define GRID_H
#include <Eigen/Dense>
#include <utility>
#include "ndim_const.h"
#include <iostream>


using namespace Eigen;
using namespace std;

class Grid
{
	Array<double, GLOBAL_NDIM, 1> center;

public:

	Array<double, GLOBAL_NDIM, 1> delta;

	Grid(Array<double, GLOBAL_NDIM, 1> in_delta, Array<double, GLOBAL_NDIM, 1> in_center);

	Array<int, GLOBAL_NDIM, 1> get_point(const Array<double, GLOBAL_NDIM, 1>& x) const;

	Array<double, Dynamic, GLOBAL_NDIM> map2x(const Array<int, Dynamic, GLOBAL_NDIM>& point);
};


#endif