#pragma once

#include "ChebTools/ChebTools.h"

// Data for a given interval in [wmin, wmax]
struct WInterval {
	const std::vector<double> wnodes;
	const std::vector<ChebTools::ChebyshevCollection> expsL, expsV; // In Theta, at each C-L node
	double wmin() const { return wnodes.back(); };
	double wmax() const { return wnodes.front(); };
};

struct CompleteWInterval {
	const std::vector<double> Wedges; 
	const std::vector<WInterval> intervals;
};

template<std::size_t N>
auto buildmat() {
    Eigen::Matrix<double, N + 1, N + 1> L; ///< Matrix of coefficients
    for (int j = 0; j <= N; ++j) {
        for (int k = j; k <= N; ++k) {
            double p_j = (j == 0 || j == N) ? 2 : 1;
            double p_k = (k == 0 || k == N) ? 2 : 1;
            L(j, k) = 2.0 / (p_j * p_k * N) * cos((j * EIGEN_PI * k) / N);
            // Exploit symmetry to fill in the symmetric elements in the matrix
            L(k, j) = L(j, k);
        }
    }
    return L;
}
const static Eigen::MatrixXd V16 = buildmat<16>();

#include "domain.cpp"