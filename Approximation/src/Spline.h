#pragma once
#ifndef spline_h
#define spline_h

#include <vector>
#include "Point.h"

namespace Com_Methods
{

	class Spline
	{
    public:
        virtual ~Spline() = default;

	public:
		virtual void Update_Spline(const std::vector<Point> &Points, const std::vector<double> &F_Value) = 0;
		virtual void Get_Value(const Point &P, double * Res)const = 0;
	};
}

#endif