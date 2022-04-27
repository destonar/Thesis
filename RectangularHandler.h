#ifndef RECTANGULARHANDLER_H
#define RECTANGULARHANDLER_H
#include <Eigen/Dense>
#include "Grid.h"
#include "ISetHandler.h"
#include "ndim_const.h"
#include <iostream>

using namespace Eigen;
using namespace std;

class RectangularHandler final: public ISetHandler
{
	Array<double, GLOBAL_NDIM, 2> bounds;
public:
	explicit RectangularHandler(const Array<double, GLOBAL_NDIM, 2>& in_bounds);

	ArrayXNi project(Grid grid) override;

	Array<double, GLOBAL_NDIM, 1> support_function(Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x) override;

	Array<double, Dynamic, 1> support_function(Array<double, Dynamic, GLOBAL_NDIM> x) override;
};
#endif