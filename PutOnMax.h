#ifndef PUTONMAX_H
#define PUTONMAX_H
#include "IOption.h"


using namespace Eigen;
using namespace std;

class PutOnMax final : public IOption
{
	double K;
public:
	explicit PutOnMax(double in_K);

	ArrayXd payoff(ArrayXNd S) override;
};


#endif