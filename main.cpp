// STD
#include <iostream>
#include <set>
#include <utility>
#include <vector>
#include <ctime>
#include <fstream>

// INTERNAL HEADERS
#include "RectangularHandler.h"
#include "Grid.h"
#include "PutOnMax.h"
#include "OptionPricer.h"
#include "NonNegativeSpaceHandler.h"
#include "EllipseHandler.h"
#include "ISetHandler.h"

// EIGEN
#include <Eigen/Dense>

#include "CallOnMax.h"
#include "CallOnMin.h"

// NAMESPACE
using namespace Eigen;
using namespace std;

//GLOBAL DIMENSION CONSTANT
#include "ndim_const.h"

// TYPEDEF FOR PAIR
typedef pair<Array<double, Dynamic, 1>*, ArrayXNi*> ValueFunc;

// MAIN FUNCTION
int main()
{
	// INPUTS
	Array<double, GLOBAL_NDIM, 2> dx;
	Array<double, GLOBAL_NDIM, 1> delta(0.08, 0.08), center(0, 0), x0(3, 3);
	Array2d mu(0, 0), sigma(0.2, 0.2);
	Matrix2d gamma;

	ArrayXXd X_01 = RowVectorXd::LinSpaced(10, 2.7, 4.5).replicate(10, 1);
	ArrayXXd X_02 = VectorXd::LinSpaced(10, 2.7, 4.5).replicate(1, 10);

	double rho = 0;

	gamma = (sigma.matrix()).asDiagonal();
	gamma(0, 1) = rho * sqrt(gamma(0, 0) * gamma(1, 1));
	gamma(1, 0) = rho * sqrt(gamma(0, 0) * gamma(1, 1));
		
	dx << -0.2, 0.2,
		-0.2, 0.2;

	//dx << 0.75, 1.25,
	//	0.75, 1.25;

	double strike = 1.5;

	//  OBJECTS OF CLASSES
	ISetHandler * supp = new RectangularHandler(dx);
	//ISetHandler * supp = new EllipseHandler(mu, gamma, 0.6);
	Grid grid(delta, center);
	IOption* opt = new CallOnMax(strike);
	ISetHandler* constraint_set = new NonNegativeSpaceHandler;

	constexpr char qh_opts = 'Qt';
	
	bool complex_initialization = false;

	if (!complex_initialization)
	{
		int N = 2;
		/*vector<int> N_iterator{1, 2, 3};

		ofstream fout;
		fout.open("output.txt");

		for (auto& k : N_iterator)
		{

			// TIMER
			clock_t timer = clock();

			// PRICER
			//OptionPricer pricer(grid, N, opt, x0, supp, constraint_set, qh_opts);
			OptionPricer pricer(grid, k, opt, x0, supp, constraint_set, qh_opts);

			// PRICER.EVALUATE
			ValueFunc res = pricer.evaluate();

			// STOP THE TIMER
			timer = clock() - timer;

			fout << "For N = " << k << " it took " << (float)timer / CLOCKS_PER_SEC << " seconds" << endl;
			cout << "For N = " << k << " it took " << (float)timer / CLOCKS_PER_SEC << " seconds" << endl;

			delete[] res.first;
			delete[] res.second;


		}*/

		// TIMER
		clock_t timer = clock();

		// PRICER
		OptionPricer pricer(grid, N, opt, x0, supp, constraint_set, qh_opts);

		// PRICER.EVALUATE
		ValueFunc res = pricer.evaluate();

		// STOP THE TIMER
		timer = clock() - timer;

		// OUTPUT
		
		cout << "V:" << endl;
		for (int t = 0; t <= N; t++)
		{
			cout << res.first[t] << endl;
		}
		cout << "___________" << endl;

		cout << "Vp:" << endl;
		for (int t = 0; t <= N; t++)
		{
			cout << res.second[t] << endl;
		}
		cout << "___________" << endl;
		
		cout << "It took " << static_cast<float>(timer) / CLOCKS_PER_SEC << " seconds" << endl;
		
		// OUTPUT IN A FILE
		/*
		ofstream fout;
		fout.open("output.txt");
		
		fout << "V:" << endl;
		for (int t = 0; t <= N; t++)
		{
			for (int i = 0; i < res.first[t].rows(); i++)
			{
				fout << res.first[t].row(i) << ", ";
			}
		}
		fout << endl << "___________" << endl;

		fout << "Vp:" << endl;

		fout << "1st column:" << endl;
		for (int t = 0; t <= N; t++)
		{
			for (int i = 0; i < res.second[t].rows(); i++)
			{
				fout << res.second[t](i, 0) << ", ";
			}
		}

		fout << endl << "2nd column:" << endl;
		for (int t = 0; t <= N; t++)
		{
			for (int i = 0; i < res.second[t].rows(); i++)
			{
				fout << res.second[t](i, 1) << ", ";
			}
		}
		fout << endl << "___________" << endl;

		fout << "It took " << (float)timer / CLOCKS_PER_SEC << " seconds" << endl;


		fout.close();
		*/
		//delete[] res.first;
		//delete[] res.second;
	}
	else
	{
		// VALUE FUNC
		ArrayXXd V(X_01.rows(), X_01.cols());
		vector<int> N_iterator{ 3};
		ofstream fout;
		fout.open("output.txt");

		// TIMER
		clock_t timer = clock();

		for (auto & k : N_iterator)
		{
			cout << "______________________" << endl;
			cout << "Calculating for N = " << k << endl;
			cout << "______________________" << endl;
			for (int i = 0; i < X_01.rows(); i++)
			{
				for (int j = 0; j < X_01.cols(); j++)
				{
					cout << "Iteration " << X_01.rows() * i + j + 1 << "/" << X_01.cols() * X_01.rows() << endl;
					// X_0 FOR EACH EVALUATION CYCLE
					Array<double, GLOBAL_NDIM, 1> X_0;
					X_0 << X_01(i, j), X_02(i, j);

					// PRICER
					OptionPricer pricer(grid, k, opt, X_0, supp, constraint_set, qh_opts);

					// PRICER.EVALUATE
					ValueFunc tmp_res = pricer.evaluate();

					// GET RESULTS
					V(i, j) = tmp_res.first[0](0, 0);

					delete[] tmp_res.first;
					delete[] tmp_res.second;
				}
			}

			// STOP THE TIMER
			timer = clock() - timer;

			// OUTPUT
			for (int i = 0; i < X_01.rows(); i++)
			{
				for (int j = 0; j < X_01.cols(); j++)
				{
					fout << V(i, j);
					if (j == X_01.cols() - 1)
						fout << ";";
					else
						fout << ", ";
				}
				fout << endl;
			}
			fout << "It took " << static_cast<float>(timer) / CLOCKS_PER_SEC << " seconds" << endl;
			
		}
		fout.close();
	}


	delete supp;
	delete constraint_set;

	// RETURN 0
	return 0;
}