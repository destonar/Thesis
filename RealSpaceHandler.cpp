#include "RealSpaceHandler.h"

Array<double, Dynamic, 1> RealSpaceHandler::support_function(Array<double, Dynamic, GLOBAL_NDIM> x)
{
	Array<double, Dynamic, 1> res(x.rows(), 1);

	for (int i = 0; i < x.rows(); i++)
	{
		ArrayXi check = (x.row(i) != 0).cast<int>();
		if (check.sum() > 0)
		{
			res(i, 0) = INFINITY;
			break;
		}
	}


	return res;
}

ArrayXNi RealSpaceHandler::project(Grid grid)
{
	throw exception("Not implemented");
}

Array<double, GLOBAL_NDIM, 1> RealSpaceHandler::support_function(Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x)
{
	throw exception("Not implemented");
}