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
	Array<double, GLOBAL_NDIM, 1> mu, EigVal;
	Array<double, GLOBAL_NDIM, GLOBAL_NDIM> sigma, U, sigma_inv;
	Matrix<double, GLOBAL_NDIM, GLOBAL_NDIM> L;


public:

	explicit EllipseHandler(const Array<double, GLOBAL_NDIM, 1>& in_mu, const Array<double, GLOBAL_NDIM, GLOBAL_NDIM>& in_sigma, double in_conf_level);

	ArrayXNi project(Grid grid) override;

	Array<double, GLOBAL_NDIM, 1> support_function(Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x) override;

	Array<double, Dynamic, 1> support_function(Array<double, Dynamic, GLOBAL_NDIM> x) override;

};
#endif

