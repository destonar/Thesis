#ifndef MINKOPERATIONS_H
#define MINKOPERATIONS_H
#include <Eigen/Dense>
#include <set>

using namespace Eigen;
using namespace std;

class MinkOperations
{
	Array<int, Dynamic, 2> x1;
	Array<int, Dynamic, 2> x2;
public:
	MinkOperations(const Array<int, Dynamic, 2>& in_x1, const Array<int, Dynamic, 2>& in_x2);

	Array<int, Dynamic, 2> minksum();
};


#endif