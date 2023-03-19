#ifndef REGULAR_PARTIONS_H
#define REGULAR_PARTIONS_H

#include <vector>

#include "Point.h"

using Point = Com_Methods::Point;

auto generateGrid(double start, double end, double step, double factor = 1) -> std::vector<Point>;

#endif // REGULAR_PARTIONS_H
