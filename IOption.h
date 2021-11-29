#ifndef IOPTION_H
#define IOPTION_H
#include <Eigen/Dense>
#include "Grid.h"
#include "ndim_const.h"

using namespace Eigen;
using namespace std;

typedef Array<double, Dynamic, GLOBAL_NDIM> ArrayXNd;

class IOption
{
public:
	virtual ~IOption() = default;
	virtual ArrayXd payoff(ArrayXNd S) = 0;
};



#endif