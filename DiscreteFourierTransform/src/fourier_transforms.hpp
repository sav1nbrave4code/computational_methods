#include <vector>
#include <complex>

#ifndef FOURIER_TRANSFORMS_HPP
#define FOURIER_TRANSFORMS_HPP

// cmm - computational methods
namespace cmm
{
    //   dft - discrete fourier transform
    auto dft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>>;

    //   idft - inverse discrete fourier transform
    auto idft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>>;

    //   fft - fast fourier transform (for even size vector)
    auto fft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>>;

    //   ifft - inverse fast fourier transform (for even size vector)
    auto ifft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>>;

    constexpr double min_border {1e-6};

} //  cmm

#endif // FOURIER_TRANSFORMS_HPP
