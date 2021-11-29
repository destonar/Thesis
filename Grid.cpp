#include"Grid.h"

Grid::Grid(Array2d in_delta, Array2d in_center)
{
	delta = std::move(in_delta);
	center = std::move(in_center);
}

Array2i Grid::get_point(const Array2d& x) const
{
	return round((x - center) / delta).cast<int>();
}

Array<double, Dynamic, 2> Grid::map2x(const Array<int, Dynamic, 2>& point)
{
	return ((point.cast<double>()).rowwise() * delta.transpose()).rowwise() - center.transpose();
}