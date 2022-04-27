#include"Grid.h"

Grid::Grid(Array<double, GLOBAL_NDIM, 1> in_delta, Array<double, GLOBAL_NDIM, 1> in_center)
{
	delta = std::move(in_delta);
	center = std::move(in_center);
}

Array<int, GLOBAL_NDIM, 1> Grid::get_point(const Array<double, GLOBAL_NDIM, 1>& x) const
{
	Array<int, GLOBAL_NDIM, 1> point = round((x - center) / delta).cast<int>();
	for (int i = 0; i < GLOBAL_NDIM; i++)
	{
		if (point(i, 0) % 2 != 0)
		{
			point(i, 0) = point(i, 0) - point(i, 0) % 2;
		}
	}

	return point;
}

Array<double, Dynamic, GLOBAL_NDIM> Grid::map2x(const Array<int, Dynamic, GLOBAL_NDIM>& point)
{
	return ((point.cast<double>()).rowwise() * delta.transpose()).rowwise() - center.transpose();
}