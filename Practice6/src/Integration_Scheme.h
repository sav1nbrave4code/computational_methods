#pragma once
#include "Point.h"
#include <vector>

namespace Com_Methods
{
	class Integration_Scheme
	{
	protected:
		std::vector<Point> Points;
		std::vector<double> Weight;
	public:
		enum Integration_Scheme_Type
		{
			Gauss1 = 1,
            Gauss2,
            Gauss3,
            Gauss4,
            Gauss5,
            Trapeze,
            Parabola
		};
	};
}

