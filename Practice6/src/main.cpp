﻿#include <iostream>
#include <cmath>

#include "Point.h"
#include "Integration_Scheme_Interval.h"

int main()
{
	//подынтегральная функция f(x) = -e^(x/2) + x
	std::function<double(const Com_Methods::Point &P)> f =
	[](const Com_Methods::Point &P) { return -std::exp(P.x() / 2) + P.x(); };
	//первообразная F(x) = x^2 / 2 - 2e^(-x/2)
	std::function<double(const Com_Methods::Point &P)> F =
	[](const Com_Methods::Point &P) { return (std::pow(P.x(), 2) / 2. -2.*std::exp(P.x()/2)); };

	//квадратурная формула Гаусс-2
	Com_Methods::Integration_Scheme_Interval Quadrature_Formula(Com_Methods::Integration_Scheme::Parabola);

	//начало и конец отрезка интегрирования
	auto Begin = Com_Methods::Point(0, 0, 0);
	auto End   = Com_Methods::Point(1, 0, 0);

	//число сегментов
	const int Num_Segments = 16;
	
	//точное значение интеграла (ф. Ньютона-Лейбница)
	double I_True = F(End) - F(Begin);

	//численное значение интеграла
	double I = Quadrature_Formula.Calculate_Integral(Begin, End, Num_Segments, f);

	std::cout << std::scientific;
	std::cout << "h = " << (End.x() - Begin.x()) / Num_Segments << std::endl;
	std::cout << "I = " << I << std::endl;
    std::cout << "I_True = " << I_True << std::endl;
	std::cout << "|I - I_True| = " << std::fabs(I - I_True) << std::endl;
}