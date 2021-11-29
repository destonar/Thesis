#include "RectangularHandler.h"

RectangularHandler::RectangularHandler(const Array22d& in_bounds)
{
	bounds = in_bounds;
}

ArrayXNi RectangularHandler::project(const Grid grid)
{
	
	Array2i r1, r2;

	r1 << grid.get_point(bounds.row(0));
	r2 << grid.get_point(bounds.row(1));

	const int ar1_size = abs(r1(0, 0) - r1(1, 0)) + 1;
	const int ar2_size = abs(r2(0, 0) - r2(1, 0)) + 1;
	Array<int, Dynamic, 2> bnd(ar1_size * ar2_size, 2);

	bnd << Matrix<int, Dynamic, 1>::LinSpaced(ar2_size * ar1_size, r1(0, 0), r1(1, 0)),
		Matrix<int, Dynamic, 1>::LinSpaced(ar2_size, r2(0, 0), r2(1, 0)).replicate(ar1_size, 1); // может быть ошибка в размерах
	return bnd;
}

Array<double, GLOBAL_NDIM, 1> RectangularHandler::support_function(const Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x)
{
	Array22d tmp = (x * bounds.col(0).transpose().replicate(2, 1)).max(x * bounds.col(1).transpose().replicate(2, 1));
	return tmp.rowwise().sum();
}

Array<double, Dynamic, 1> RectangularHandler::support_function(Array<double, Dynamic, GLOBAL_NDIM> x)
{
	throw exception("Not implemented");
}

