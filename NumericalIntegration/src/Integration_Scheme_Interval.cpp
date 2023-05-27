#include <cmath>
#include <utility>

#include "Integration_Scheme_Interval.h"

namespace Com_Methods
{
	Integration_Scheme_Interval::Integration_Scheme_Interval(Integration_Scheme_Type Type)
	{
		switch (Type)
		{
            case Gauss2:
            {
                Weight = {1, 1};
                Points = {Point{-1. / std::sqrt(3), 0, 0}, Point{1. / std::sqrt(3), 0, 0}};
                break;
            }
            case Parabola:
            {
                Weight = { 1./3, 4./3, 1./3};
                Points = {Point{-1,0 ,0}, Point{0, 0, 0}, Point{1, 0, 0}};
                break;
            }
            default:
            {
                throw "Unknown scheme type";
            }
		}
    }


	double Integration_Scheme_Interval:: Calculate_Integral(
								         const Point &Begin,
								         const Point &End,
								         int Number_Segments,
								         const std::function<double(const Point &P)>&Func) const
	{
		double Result = 0.0;
		double X0;
		double h = (End.x() - Begin.x()) / Number_Segments;
		for (int i = 0; i < Number_Segments; i++)
		{
			X0 = Begin.x() + i * h;
			for (int Integ_Point = 0; Integ_Point < Points.size(); Integ_Point++)
			{
				auto P = Point(X0 + (1 + Points[Integ_Point].x()) * h / 2.0, 0, 0);
				Result += Weight[Integ_Point] * Func(P);
			}
		}
		return Result * (h / 2.0);
	}
}
