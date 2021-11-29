#include "EllipseHandler.h"


EllipseHandler::EllipseHandler(const Array<double, 2, 1>& in_mu, const Array<double, 2, 2>& in_sigma,
	const double in_conf_level)
{
	mu = in_mu;
	sigma = in_sigma * (-2 * log(1 - in_conf_level));
	const EigenSolver<Matrix2d> es(sigma.matrix(), true);
	EigVal = es.eigenvalues().real();
	U = es.eigenvectors().real();
	L = (EigVal.matrix()).asDiagonal();
	L = sqrt(L.array());
	sigma_inv = sigma.matrix().fullPivLu().solve(Matrix<double, 2, 2>::Identity());
}

ArrayXNi EllipseHandler::project(Grid grid)
{
	const double R = L.diagonal().cwiseAbs().maxCoeff() + grid.delta.maxCoeff();
	Array<double, GLOBAL_NDIM, 2> bounds;
	bounds << mu - R, mu + R;

	ISetHandler* tmp_supp = new RectangularHandler(bounds);

	ArrayXNi S = tmp_supp->project(grid);

	vector<int> ind;

	ArrayXNd tmpSmu = grid.map2x(S) - mu.transpose().replicate(S.rows(), 1);

	ArrayXd r2 = ((tmpSmu.matrix() * sigma_inv.matrix()).array() * tmpSmu).rowwise().sum();

	for (int i = 0; i < r2.rows(); i++)
	{
		if (r2(i, 0) <= 1)
			ind.push_back(i);
	}

	delete tmp_supp;

	return S(ind, all);
}

Array<double, GLOBAL_NDIM, 1> EllipseHandler::support_function(Array<double, GLOBAL_NDIM, GLOBAL_NDIM> x)
{
	Array22d tmp = L.matrix() * U.transpose().matrix() * x.transpose().matrix();
	return ((abs(tmp).pow(2)).colwise().sum()).pow(0.5);
	//return tmp.matrix().colwise().norm();
}

Array<double, Dynamic, 1> EllipseHandler::support_function(Array<double, Dynamic, GLOBAL_NDIM> x)
{
	throw exception("Not implemented");
}
