#pragma once
#ifndef point_h
#define point_h

namespace Com_Methods
{

	class Point
	{
	private:

		double X, Y, Z;
	public:
		explicit Point(double x = 0, double y = 0, double z = 0);

		[[nodiscard]] double x() const;
		[[nodiscard]] double y() const;
		[[nodiscard]] double z() const;
	};
}

#endif

