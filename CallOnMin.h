#ifndef CALLONMIN_H
#define CALLONMIN_H
#include "IOption.h"


using namespace Eigen;
using namespace std;

class CallOnMin final : public IOption
{
	double K;
public:
	explicit CallOnMin(double in_K);

	ArrayXd payoff(ArrayXNd S) override;
};


#endif