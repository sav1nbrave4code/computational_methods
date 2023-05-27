#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Point.h"
#include "Spline.h"
#include "Cubic_Interpolation_Spline_1D.h"
#include "regular_partions.h"

using Interface = Com_Methods::Spline;
using Cubic = Com_Methods::Cubic_Interpolation_Spline_1D;
using Point = Com_Methods::Point;

auto main() -> int
{
    std::ofstream              out    {"out.txt"};
    std::unique_ptr<Interface> spline {new Cubic};
    std::unique_ptr<double>    result {new double[3]{}};

    std::cout.rdbuf(out.rdbuf());

    auto function {[](const double x) -> double{
        return 1 / x + 0.1;
    }};
    auto functionFirstDerivative {[](const double x){
        return -1 / (x * x);
    }};
    auto functionSecondDerivative {[](const double x){
        return 2 / (x * x * x);
    }};

    double              step     {0.25};
    std::vector<Point>  grid     {generateGrid(0.5, 1, step)};
    std::vector<Point>  testGrid {Point{0.51}, Point{0.57}, Point{0.62}, Point{0.69}, Point{0.73},
                                  Point{0.78}, Point{0.84}, Point{0.87}, Point{0.92}, Point{0.98}};
    std::vector<double> values   {};

    std::cout << std::setprecision(6);

    std::cout << "function analyses x, f(x), f'(x), f''(x)" << std::endl;
    for (const auto& elem : testGrid)
    {
        std::cout << elem.x() << " " << function(elem.x()) << " " << functionFirstDerivative(elem.x()) << " "
        << functionSecondDerivative(elem.x()) << std::endl;
    }
    std::cout << std::endl;

    for (const auto& elem : grid)
    {
        values.push_back(function(elem.x()));
    }
    spline->Update_Spline(grid, values);
    std::cout << "spline analyses h=" << step << std::endl;
    std::cout << "x\t\ts(xi)\t\ts'(xi)\t\ts''(xi)" << std::endl;
    for (auto i {0}; i < grid.size(); ++i)
    {
        spline->Get_Value(grid[i], result.get());

        std::cout << grid[i].x() << " " << result.get()[0] << " " << result.get()[1] << " " << result.get()[2];
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "testing on testGrid s(x), s'(x)" << std::endl;
    for (const auto& elem : testGrid)
    {
        spline->Get_Value(elem, result.get());
        std::cout << elem.x() << " " << result.get()[0] << " " << result.get()[1] << " " << result.get()[2];
        std::cout << std::endl;
    }

    grid = generateGrid(0.5, 1, step / 2.);
    values.clear();
    for (const auto& elem : grid)
    {
        values.push_back(function(elem.x()));
    }
    spline->Update_Spline(grid, values);
    std::cout << "spline analyses h/2=" << step / 2. << std::endl;
    std::cout << "x\t\ts(xi)\t\ts'(xi)\t\ts''(xi)" << std::endl;
    for (auto i {0}; i < grid.size(); ++i)
    {
        spline->Get_Value(grid[i], result.get());

        std::cout << grid[i].x() << " " << result.get()[0] << " " << result.get()[1] << " " << result.get()[2];
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "testing on testGrid s(x), s'(x)" << std::endl;
    for (const auto& elem : testGrid)
    {
        spline->Get_Value(elem, result.get());
        std::cout << elem.x() << " " << result.get()[0] << " " << result.get()[1] << " " << result.get()[2];
        std::cout << std::endl;
    }

    grid = generateGrid(0.5, 1, step / 4.);
    values.clear();
    for (const auto& elem : grid)
    {
        values.push_back(function(elem.x()));
    }
    spline->Update_Spline(grid, values);
    std::cout << "spline analyses h/2=" << step / 4. << std::endl;
    std::cout << "x\t\ts(xi)\t\ts'(xi)\t\ts''(xi)" << std::endl;
    for (auto i {0}; i < grid.size(); ++i)
    {
        spline->Get_Value(grid[i], result.get());

        std::cout << grid[i].x() << " " << result.get()[0] << " " << result.get()[1] << " " << result.get()[2];
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "testing on testGrid s(x), s'(x)" << std::endl;
    for (const auto& elem : testGrid)
    {
        spline->Get_Value(elem, result.get());
        std::cout << elem.x() << " " << result.get()[0] << " " << result.get()[1] << " " << result.get()[2];
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "spline values at critical points" << std::endl;
    spline->Get_Value(Point{0.5}, result.get());
    std::cout << "x: 0.5 " << "f(x): " << result.get()[0] << " f'(x): " << result.get()[1] << std::endl;
    spline->Get_Value(Point{1}, result.get());
    std::cout << "x: 1 " << "f(x): " << result.get()[0] << " f'(x): " << result.get()[1] << std::endl;
}