#include"PutOnMax.h"

PutOnMax::PutOnMax(const double in_K)
{
	K = in_K;
}

ArrayXd PutOnMax::payoff(ArrayXNd S)
{
	const ArrayXd tmp = K - S.rowwise().maxCoeff();
	return tmp.max(ArrayXd::Zero(tmp.rows(), 1));
}
