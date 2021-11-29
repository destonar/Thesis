#ifndef GRID_H
#define GRID_H
#include <Eigen/Dense>
#include <utility>


using namespace Eigen;
using namespace std;

class Grid
{
	Array2d center;

public:

	Array2d delta;

	Grid(Array2d in_delta, Array2d in_center);

	Array2i get_point(const Array2d& x) const;

	Array<double, Dynamic, 2> map2x(const Array<int, Dynamic, 2>& point);
};


#endif