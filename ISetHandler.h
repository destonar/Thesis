#ifndef ISETHANDLER_H
#define ISETHANDLER_H
#include <Eigen/Dense>
#include "Grid.h"
#include "ndim_const.h"

using namespace Eigen;
using namespace std;

typedef Array<double, Dynamic, GLOBAL_NDIM> ArrayXNd;
typedef Array<int, Dynamic, GLOBAL_NDIM> ArrayXNi;

class ISetHandler
{
public:
	virtual ~ISetHandler() = default;
	virtual ArrayXNi project(Grid grid) = 0;
	virtual Array<double, GLOBAL_NDIM, 1> support_function(Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x) = 0;
	virtual Array<double, Dynamic, 1> support_function(Array<double, Dynamic, GLOBAL_NDIM> x) = 0;
};



#endif