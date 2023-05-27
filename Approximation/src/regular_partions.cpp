#include <exception>

#include "regular_partions.h"

auto generateGrid(const double start, const double end, double step, const double factor) -> std::vector<Point>
{
    if (step <= 0. || start >= end || factor < 1.)
    {
        throw std::exception{};
    }

    std::vector<Point> result {};

    for (double i {start}; i <= end + 1e-7; i += step, step *= factor)
    {
        result.emplace_back(i);
    }

    return result;
}

