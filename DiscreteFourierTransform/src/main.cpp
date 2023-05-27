#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include <numbers>
#include <chrono>
#include <fstream>

#include "fourier_transforms.hpp"

auto print_table(auto& z, auto& dft, auto N)-> void
{
    auto width {22};

    std::cout << std::left << std::setw(width) << "Number"
              << std::left << std::setw(width) << "Re(Z)"
              << std::left << std::setw(width) << "Im(Z)"
              << std::left << std::setw(width) << "Re(DFT)"
              << std::left << std::setw(width) << "Im(DFT)"
              << std::left << std::setw(width) << "Abs(DFT)"
              << std::left << std::setw(width) << "Arg(DFT)" << std::endl;

    for (auto i {0}; i < N; i++)
    {
        auto abs {std::abs(dft[i])};
        auto arg {std::arg(dft[i])};

        if (abs < cmm::min_border)
        {
            continue;
        }

        std::cout << std::left << std::setw(width) << i
                  << std::left << std::setw(width) << z[i].real()
                  << std::left << std::setw(width) << z[i].imag()
                  << std::left << std::setw(width) << dft[i].real()
                  << std::left << std::setw(width) << dft[i].imag()
                  << std::left << std::setw(width) << abs
                  << std::left << std::setw(width) << arg << std::endl;
    }
}


template<class T>
auto write_values_to_file(std::vector<T>& values, std::ofstream& file) -> void;

template<>
auto write_values_to_file<std::complex<double>>(std::vector<std::complex<double>>& values, std::ofstream & file) -> void
{
    for (const auto item : values)
    {
        file << std::fixed << std::setprecision(6) << item.real() << std::endl;
    }
    file.close();
}

template<>
auto write_values_to_file<double>(std::vector<double>& values, std::ofstream & file) -> void
{
    for (const auto item : values)
    {
        file << std::fixed << std::setprecision(6) << item << std::endl;
    }
    file.close();
}

auto main() -> int
{
	int N = 512; // vector space size
	std::vector<std::complex<double>> Z(N); // data

    //  19th variant data
    double A   {-40};
    double B   {0.76};
    double w   {135};
    double phi {std::numbers::pi / 4};

    //  fill data
	for (int i = 0; i < N; i++)
	{
        Z[i] = { A + B * std::cos(2 * std::numbers::pi * w * i / N + phi), 0};
	}

    auto begin {std::chrono::high_resolution_clock::now()};
	auto dft_result {cmm::dft(Z)};
    auto end {std::chrono::high_resolution_clock::now()};
    std::cout << "DFT compute time: " << std::chrono::duration<double>(end - begin).count() << std::endl;

    begin = std::chrono::high_resolution_clock::now();
    auto fft_result {cmm::fft(Z)};
    end = std::chrono::high_resolution_clock::now();
    std::cout << "FFT compute time: " << std::chrono::duration<double>(end - begin).count() << std::endl;
    std::cout << std::endl;


    std::cout << "Original signal" << std::endl;
    print_table(Z, dft_result, N);
    std::cout << std::endl;

    std::ofstream non_filtered_file {"non-filtered.txt"};
    std::ofstream filtered_file {"filtered.txt"};
    std::ofstream range_file {"range.txt"};

    std::vector<double> range {};
    range.resize(N);

    for (auto i {0}; i < N; ++i)
    {
        range[i] = 2 * std::numbers::pi / N * i;
    }

    //  non-filtered signal
    for (auto i {0}; i < N; ++i)
    {
        Z[i] = {std::cos(2 * std::numbers::pi * i / N) + 0.01 * std::cos( 2 * std::numbers::pi * i * w / N), 0};
    }
    std::cout << "Non-filtered signal" << std::endl;
    dft_result = cmm::fft(Z);
    print_table(Z, dft_result, N);

    //  garbage waves
    dft_result[135] = 0;
    dft_result[377] = 0;

    auto idft_values {cmm::idft(dft_result)};

    write_values_to_file<double>(range, range_file);
    write_values_to_file<std::complex<double>>(Z, non_filtered_file);
    write_values_to_file<std::complex<double>>(idft_values, filtered_file);
}