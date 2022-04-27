#include "MinkOperations.h"

MinkOperations::MinkOperations(const Array<int, Dynamic, GLOBAL_NDIM>& in_x1, const Array<int, Dynamic, GLOBAL_NDIM>& in_x2)
{
	x1 = in_x1;
	x2 = in_x2;
}

Array<int, Dynamic, GLOBAL_NDIM> MinkOperations::minksum()
{
	set<vector<int>> s;

	for (int i = 0; i < x1.rows(); i++)
	{
		for (int j = 0; j < x2.rows(); j++)
		{
			vector<int> vec_row_sum;
			vec_row_sum.reserve(GLOBAL_NDIM);
			for (int k = 0; k < GLOBAL_NDIM; k++)
			{
				vec_row_sum.push_back(x1(i, k) + x2(j, k));
			}

			if (s.find(vec_row_sum) == s.end())
			{
				s.insert(vec_row_sum);
			}

		}
	}
	Array<int, Dynamic, GLOBAL_NDIM> sum(s.size(), GLOBAL_NDIM);

	int i = 0;
	for (auto& vec : s)
	{
		for (int j = 0; j < GLOBAL_NDIM; j++)
		{
			sum(i, j) = vec[j];
		}
		i++;
	}
	return sum;
}
