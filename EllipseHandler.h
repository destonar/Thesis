#ifndef ELLIPSEHANDLER_H
#define ELLIPSEHANDLER_H
#include "Grid.h"
#include "ISetHandler.h"
#include "RectangularHandler.h"
#include <vector>

using namespace Eigen;
using namespace std;


class EllipseHandler final : public ISetHandler
{
	Array2d mu, EigVal;
	Array22d sigma, U, sigma_inv;
	Matrix2d L;


public:

	explicit EllipseHandler(const Array2d& in_mu, const Array22d& in_sigma, double in_conf_level);

	ArrayXNi project(Grid grid) override;

	Array<double, GLOBAL_NDIM, 1> support_function(Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x) override;

	Array<double, Dynamic, 1> support_function(Array<double, Dynamic, GLOBAL_NDIM> x) override;

};
#endif


