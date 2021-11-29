#ifndef OPTIONPRICER_H
#define OPTIONPRICER_H
// STD LIBS
#include <iostream>

// EIGEN
#include <Eigen/Dense>

// CPPIPM
#include "cppipm.h"
#include "mpsReader.h"
#include "include_libs.h"

// QHULL
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullUser.h"

// INTERNAL HEADERS
#include "Grid.h"
#include "ISetHandler.h"
#include "PutOnMax.h"
#include "MinkOperations.h"
#include "ndim_const.h"

// DEBUG OUTPUT
#include <fstream>






using AlgorithmPtr = std::unique_ptr<Algorithm>;

typedef Array<double, Dynamic, GLOBAL_NDIM> ArrayXNd;
typedef Array<int, Dynamic, GLOBAL_NDIM> ArrayXNi;
typedef pair<Array<double, Dynamic, 1>*, ArrayXNi*> ValueFunc;

class OptionPricer
{
	Grid grid;
	int N;
	PutOnMax option;
	Array2d x0;
	ISetHandler* price_support;
	//EllipseHandler price_support;
	Array2i p0;
	ISetHandler* constraint_set;
	const char qh_opts;
	ofstream fout;



public:

	OptionPricer(Grid in_grid, int in_N, PutOnMax in_option, Array2d in_x0, ISetHandler* in_price_support, ISetHandler* in_constraint_set, const char in_qh_opts);

	ArrayXNi get_support_set(const ArrayXNi& curr_set, ISetHandler* increment) const;

	ArrayXNi* generate_evaluation_point_lists(Array2i p_0, ISetHandler* pricer_support, int N) const;

	VectorXd get_max_coordinates(ArrayXNd x, const ArrayXd& f, Array<double, 1, GLOBAL_NDIM> z) const;

	double find_u(ArrayXNd x, ArrayXd v, Array<double, 1, GLOBAL_NDIM> z) const;

	double find_rho(ArrayXNd x, ArrayXd v, ArrayXNd K_x, ArrayXNd convdK_x) const;

	ValueFunc evaluate();
};



#endif