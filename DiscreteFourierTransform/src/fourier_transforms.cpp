#include "fourier_transforms.hpp"

namespace cmm {

    auto dft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>>
    {
        auto N {data.size()};
        std::complex<double> Exp {};
        std::vector<std::complex<double>> result {};

        result.resize(N);

        for (int n = 0; n < N; n++)
        {
            for (int m = 0; m < N; m++)
            {
                Exp = std::complex<double>{std::cos(2 * std::numbers::pi * m * n / static_cast<double>(N)),
                                           std::sin(2 * std::numbers::pi * m * n / static_cast<double>(N))};
                result[m] += data[n] * Exp;
            }
        }

        return result;
    }

    auto idft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>>
    {
        auto N {data.size()};
        std::complex<double> Exp;
        std::vector<std::complex<double>> result {};

        result.resize(N);

        for (int n = 0; n < N; n++)
        {
            for (int m = 0; m < N; m++)
            {
                Exp = std::complex<double>{std::cos(2 * std::numbers::pi * m * n / static_cast<double>(N)),
                                           std::sin(2 * std::numbers::pi * m * n / static_cast<double>(N))};
                result[n] += data[m] * Exp;
            }
            result[n] *= (1 / double(N));
        }

        return result;
    }

    auto fft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>> {
        auto N{data.size()};
        auto M{N / 2};
        std::complex<double> Exp {};
        std::complex<double> U {};
        std::complex<double> V {};
        std::vector<std::complex<double>> result {};

        result.resize(N);

        for (int m = 0; m < M; m++) {
            U = std::complex<double>{0, 0};
            V = std::complex<double>{0, 0};
            for (auto n{0}; n < M; n++) {
                Exp = std::complex<double>{std::cos(-2.0 * std::numbers::pi * m * n / static_cast<double>(M)),
                                           std::sin(-2.0 * std::numbers::pi * m * n / static_cast<double>(M))};
                U += data[2 * n] * Exp;
                V += data[2 * n + 1] * Exp;
            }
            Exp = std::complex<double>{std::cos(-2.0 * std::numbers::pi * m / static_cast<double>(N)),
                                       std::sin(-2.0 * std::numbers::pi * m / static_cast<double>(N))};
            result[m] = U + Exp * V;
            result[m + M] = U - Exp * V;
        }

        return result;
    }

    auto ifft(const std::vector<std::complex<double>> &data) -> std::vector<std::complex<double>> {
        auto N{data.size()};
        std::complex<double> val {};
        std::vector<std::complex<double>> result {};

        result.resize(N);

        for (int i = 1; i <= N / 2; i++) {
            val = data[i];
            result[i] = data[N - i] / static_cast<double>(N);
            result[N - i] = val / static_cast<double>(N);
        }
        result[0] /= static_cast<double>(N);

        return result;
    }

} // cmm