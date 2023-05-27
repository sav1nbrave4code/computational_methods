#pragma once
#ifndef Smoothing_Spline_1D_h
#define Smoothing_Spline_1D_h

#include "Spline.h"

namespace Com_Methods
{

	class Smoothing_Spline_1D : public Spline
	{
	private:

		double SMOOTH;

		std::vector<Point> Points;

		std::vector<double> alpha;


		void Transition_To_Master_Element(int Seg_Num, const double &X, double &Ksi) const;


		double Basis_Function(int Number, const double &Ksi) const;


		double Der_Basis_Function(int Number, const double & Ksi) const;
	public:

		Smoothing_Spline_1D(const double &SMOOTH);

		void Update_Spline(const std::vector<Point> &Points, const std::vector<double> &F_Value) override;

		void Get_Value(const Point &P, double * Res)const override;
	};
}

#endif