#include"CallOnMax.h"

CallOnMax::CallOnMax(const double in_K)
{
	K = in_K;
}

ArrayXd CallOnMax::payoff(ArrayXNd S)
{
	const ArrayXd tmp = S.rowwise().maxCoeff() - K;
	return tmp.max(ArrayXd::Zero(tmp.rows(), 1));
}
