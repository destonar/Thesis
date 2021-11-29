#ifndef CALLONMAX_H
#define CALLONMAX_H
#include "IOption.h"


using namespace Eigen;
using namespace std;

class CallOnMax final : public IOption
{
	double K;
public:
	explicit CallOnMax(double in_K);

	ArrayXd payoff(ArrayXNd S) override;
};


#endif