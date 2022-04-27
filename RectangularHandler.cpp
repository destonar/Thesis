#include "RectangularHandler.h"

RectangularHandler::RectangularHandler(const Array<double, GLOBAL_NDIM, 2>& in_bounds)
{
	bounds = in_bounds;
}

ArrayXNi RectangularHandler::project(const Grid grid)
{

	Array<int, GLOBAL_NDIM, 2> projected_bounds;

	projected_bounds << grid.get_point(bounds.transpose().row(0)), grid.get_point(bounds.transpose().row(1));

	Array<int, GLOBAL_NDIM, 1> array_sizes;
	int projection_rows = 1;

	for (int i = 0; i < GLOBAL_NDIM; i++)
	{
		const int cur_size = abs(projected_bounds(i, 0) - projected_bounds(i, 1) + 1);
		array_sizes(i) = cur_size;
		projection_rows *= cur_size;
	}

	Array<int, GLOBAL_NDIM, Dynamic> projection(GLOBAL_NDIM, projection_rows);
	int cur_repetitions = 1;
	for (int i = 0; i < GLOBAL_NDIM; i++)
	{
		if (i == 0)
		{
			projection.row(i) = Matrix<int, Dynamic, 1>::LinSpaced(projection_rows, projected_bounds(i, 0), projected_bounds(i, 1)).transpose();
		}
		else
		{
			projection_rows /= array_sizes(i - 1);
			cur_repetitions *= array_sizes(i - 1);
			projection.row(i) = Matrix<int, Dynamic, 1>::LinSpaced(projection_rows, projected_bounds(i, 0), projected_bounds(i, 1)).replicate(cur_repetitions, 1).transpose();
		}
	}

	return projection.transpose();
}

Array<double, GLOBAL_NDIM, 1> RectangularHandler::support_function(const Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x)
{
	Array<double, GLOBAL_NDIM, GLOBAL_NDIM> tmp = (x * bounds.col(0).transpose().replicate(2, 1)).max(x * bounds.col(1).transpose().replicate(2, 1));
	return tmp.rowwise().sum();
}

Array<double, Dynamic, 1> RectangularHandler::support_function(Array<double, Dynamic, GLOBAL_NDIM> x)
{
	throw exception("Not implemented");
}


