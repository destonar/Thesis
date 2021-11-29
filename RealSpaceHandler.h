#ifndef REALSPACEHANDLER_H
#define REALSPACEHANDLER_H
#include <Eigen/Dense>
#include "ISetHandler.h"

using namespace Eigen;
using namespace std;

class RealSpaceHandler : public ISetHandler
{
public:

	Array<double, Dynamic, 1> support_function(Array<double, Dynamic, GLOBAL_NDIM> x) override;

	Array<double, GLOBAL_NDIM, 1> support_function(Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x) override;

	ArrayXNi project(Grid grid) override;
};

#endif
