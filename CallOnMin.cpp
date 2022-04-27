#include"CallOnMin.h"

CallOnMin::CallOnMin(const double in_K)
{
	K = in_K;
}

ArrayXd CallOnMin::payoff(ArrayXNd S)
{
	const ArrayXd tmp = S.rowwise().minCoeff() - K;
	return tmp.max(ArrayXd::Zero(tmp.rows(), 1));
}
