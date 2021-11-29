	#include "MinkOperations.h"

MinkOperations::MinkOperations(const Array<int, Dynamic, 2>& in_x1, const Array<int, Dynamic, 2>& in_x2)
{
	x1 = in_x1;
	x2 = in_x2;
}

Array<int, Dynamic, 2> MinkOperations::minksum()
{
	set<pair<int, int>> s;

	for (int i = 0; i < x1.rows(); i++)
	{
		for (int j = 0; j < x2.rows(); j++)
		{
			pair<int, int> p(x1(i, 0) + x2(j, 0), x1(i, 1) + x2(j, 1));

			if (s.find(p) == s.end())
			{
				s.insert(p);
			}

		}
	}
	Array<int, Dynamic, 2> sum(s.size(), 2);
	int i = 0;
	for (auto& p : s)
	{

		sum.row(i) << p.first, p.second;
		i++;
	}
	return sum;
}
