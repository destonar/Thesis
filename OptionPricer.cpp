#include "OptionPricer.h"

#include <utility>

using namespace Eigen;
using namespace std;
using namespace orgQhull;

OptionPricer::OptionPricer(Grid in_grid, const int in_N, PutOnMax in_option, Array2d in_x0, ISetHandler* in_price_support, ISetHandler* in_constraint_set, const char in_qh_opts) :
	grid(move(in_grid)), N(in_N), option(move(in_option)), x0(move(in_x0)), price_support(in_price_support), constraint_set(in_constraint_set), qh_opts(in_qh_opts)
{
	p0 = grid.get_point(x0);
	//fout.open("debug_output.txt");
}

ArrayXNi OptionPricer::get_support_set(const ArrayXNi& curr_set, ISetHandler* increment) const
{
	MinkOperations sets(curr_set, increment->project(grid));
	return sets.minksum();
}

ArrayXNi* OptionPricer::generate_evaluation_point_lists(Array2i p_0, ISetHandler* pricer_support, int N) const
{
	auto Vp = new ArrayXNi[1 + N];
	Vp[0] = p_0.transpose();
	for (int i = 0; i < N; i++)
	{
		Vp[1 + i] = get_support_set(Vp[i], pricer_support);
		
	}
	return Vp;
}

VectorXd OptionPricer::get_max_coordinates(ArrayXNd x, const ArrayXd& f, Array<double, 1, GLOBAL_NDIM> z) const
{
	Array<double, GLOBAL_NDIM + 1, Dynamic> A(GLOBAL_NDIM + 1, f.rows());
	ArrayXd tmp_ones = ArrayXd::Ones(f.rows(), 1);
	A << tmp_ones.transpose(), x.transpose();

	Array<double, GLOBAL_NDIM + 1, 1> b;
	b << 1, z.transpose();

	ArrayXd c = -f;

	AlgorithmPtr linprog_problem = std::make_unique<cppipm>(A, b, c);
	linprog_problem->solve();

	VectorXd res = linprog_problem->get_res();
	linprog_problem.reset(nullptr);
	return res;
}

double OptionPricer::find_u(ArrayXNd x, ArrayXd v, Array<double, 1, GLOBAL_NDIM> z) const
{
	//double const eps = numeric_limits< double >::epsilon();
	if (v.maxCoeff() - v.minCoeff() < 0.0000001)
	{
		Array<double, Dynamic, 1> tmp_max = (x - z.replicate(x.rows(), 1)).rowwise().maxCoeff();
		Array<double, Dynamic, 1>::Index minRow, minCol;
		tmp_max.maxCoeff(&minRow, &minCol);
		return v(minRow, minCol);
	}
	vector<int> resize_index_vector;
	for (int resize_index = 0; resize_index < v.rows(); resize_index++)
	{
		if (v(resize_index, 0) > 0)
		{
			resize_index_vector.push_back(resize_index);
		}
	}
	long long resize_vector_size = resize_index_vector.size();

	Array<double, Dynamic, GLOBAL_NDIM + 1, RowMajor> points(v.rows() + resize_vector_size, GLOBAL_NDIM + 1);

	ArrayXNd x_resized(x.rows() + resize_vector_size, GLOBAL_NDIM);

	Array<double, Dynamic, 1> v_resized(v.rows() + resize_vector_size, 1);
	ArrayXd tmp_zeros = ArrayXd::Zero(resize_vector_size, 1);

	x_resized << x, x(resize_index_vector, all);
	v_resized << v, tmp_zeros;

	points << x_resized, v_resized;

	//fout << "Points before qhull:" << endl << points << endl;

	Qhull qhull;
	auto qhull_points = new PointCoordinates(GLOBAL_NDIM + 1, "");

	// convert eigen array to vector<double>
	vector<double> vec_points;
	vec_points.reserve(points.size());

	for (int j = 0; j < points.size(); j++)
	{
		vec_points.push_back(*(points.data() + j));
	}
	// append

	qhull_points->append(vec_points);

	// run the algorithm
	qhull.runQhull(qhull_points->comment().c_str(), qhull_points->dimension(), qhull_points->count(),
		&*qhull_points->coordinates(), "QJ");


	// find facet vertices
	QhullFacetList facets = qhull.facetList();
	vector<int> vertices_ind;

	QhullFacetListIterator j(facets);
	while (j.hasNext()) {
		QhullFacet f = j.next();

		if (!f.isGood()) {
			// ignore facet
		}
		else {

			QhullVertexSetIterator k(f.vertices());
			while (k.hasNext()) {
				QhullVertex vertex = k.next();
				QhullPoint p = vertex.point();
				if (p.id() < v.rows())
					vertices_ind.push_back(p.id());
			}

		}
	}

	// get vector of indices
	sort(vertices_ind.begin(), vertices_ind.end());
	vertices_ind.erase(unique(vertices_ind.begin(), vertices_ind.end()), vertices_ind.end());

	// convert to eigen array
	Array<double, Dynamic, 3> res_vert = points(vertices_ind, all);

	//fout << "Points after qhull:" << endl << res_vert << endl;

	delete qhull_points;
	// get_max_coordinates
	return (get_max_coordinates(res_vert.block(0, 0, res_vert.rows(), GLOBAL_NDIM), res_vert.block(0, GLOBAL_NDIM, res_vert.rows(), 1), z)).dot(v(vertices_ind, all).matrix());

}

double OptionPricer::find_rho(ArrayXNd x, ArrayXd v, ArrayXNd K_x_in, ArrayXNd convdK_x_in) const
{

	Array<double, Dynamic, 1> csupp = constraint_set->support_function(convdK_x_in);
	vector<int> ind2;
	Array<bool, Dynamic, 1> check = csupp.isFinite();
	bool inf_flag = true;
	for (int check_for_inf_index = 0; check_for_inf_index < check.rows(); check_for_inf_index++)
	{
		if (check(check_for_inf_index) == true)
		{
			ind2.push_back(check_for_inf_index);
			inf_flag = false;
		}
	}
	if (inf_flag)
	{
		cout << "Support function is infinite" << endl;
		return INFINITY;
	}

	ArrayXNd K_x = K_x_in(ind2, all), convdK_x = convdK_x_in(ind2, all);

	const ArrayXd supp_func = csupp(ind2, all);


	ArrayXd res_u(K_x.rows(), 1);

	for (int u_index = 0; u_index < K_x.rows(); u_index++)
	{
		// step one: find u
		res_u(u_index, 0) = find_u(x, v, K_x.row(u_index));
	}

	
	ArrayXd::Index max_ind;

	return (res_u - supp_func).maxCoeff(&max_ind);
}

ValueFunc OptionPricer::evaluate()
{
	auto Vp = generate_evaluation_point_lists(p0, price_support, N);

	ArrayXNd x = grid.map2x(Vp[N]);

	auto Vf = new Array<double, Dynamic, 1>[1 + N];

	Vf[N] = option.payoff(x);


	for (int t = N - 1; t >= 0; t--) // reversed time
	{
		ArrayXd res(Vp[t].rows(), 1);


		for (int i = 0; i < Vp[t].rows(); i++) // loop by pregenerated points for V
		{
			vector<int> ind;
			
			ArrayXNi K = get_support_set(Vp[t].row(i), price_support);
			
			for (int k = 0; k < Vp[t + 1].rows(); k++) // find which points belong to K_t
				for (int m = 0; m < K.rows(); m++)
				{
					if ((Vp[t + 1].row(k) - K.row(m)).abs().sum() == 0)
					{
						ind.push_back(k);
						break;
					}
				}
			// step two: find rho
			res(i, 0) = find_rho(grid.map2x(Vp[1 + t](ind, all)), Vf[1 + t](ind, all),
				grid.map2x(K), grid.map2x(K.rowwise() - Vp[t].row(i)));
			
		}
		Vf[t] = res;
		
	}

	ValueFunc value_func(Vf, Vp);

	return value_func;
}



