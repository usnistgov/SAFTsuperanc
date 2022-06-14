#pragma once

#include "PCSAFTsuperancillary_data.hpp"
#include "PCSAFTsuperancillary_crit.hpp"

using namespace ChebTools;

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

const auto& get_interval(double w) {
    using namespace PCSAFTSuperAncillary;
    /// Return the index of the expansion that is desired
    auto get_index = [&](double w) {
        auto midpoint_Knuth = [](int x, int y) { 
            return (x & y) + ((x ^ y) >> 1); 
        };
        int iL = 0, iR = static_cast<int>(domain.Wedges.size()) - 2, iM;
        while (iR - iL > 1) {
            iM = midpoint_Knuth(iL, iR);
            if (w >= domain.intervals[iM].wmin()) {
                iL = iM;
            }
            else {
                iR = iM;
            }
        }
        return (w < domain.intervals[iL].wmax()) ? iL : iR;
    };
    int i = get_index(w);
    return domain.intervals[i];
}

auto get_funcvals(double Theta, const WInterval& interval) {
    Eigen::Index Nm = interval.expsL.size()-1;
    Eigen::ArrayXd rhoLfvals(Nm+1), rhoVfvals(Nm+1);
    for (auto i = 0; i <= Nm; ++i) {
        rhoLfvals[i] = interval.expsL[i](Theta);
        rhoVfvals[i] = interval.expsV[i](Theta);
    }
    return std::make_tuple(rhoLfvals, rhoVfvals);
}

auto get_expansion(const Eigen::ArrayXd& vals, double wmin, double wmax) {
    return ChebyshevExpansion((V16 * vals.matrix()).eval(), wmin, wmax);
}

auto PCSAFTsuperanc_rhoLV(double Ttilde, double m){
    // 
    using namespace PCSAFTSuperAncillary;
    auto Ttilde_crit = cc_Ttilde(1/m);
    auto Ttilde_min = exp(-2.20078778)*pow(m, 0.37627892)*Ttilde_crit;
    auto Theta = (Ttilde - Ttilde_min) / (Ttilde_crit - Ttilde_min);
    // Bisection to find the right interval in w=1/m
    const auto& interval = get_interval(1/m); 
    double wmin = interval.wmin();
    double wmax = interval.wmax();
    // Then evaluate the values at the nodes of w for each phase in this interval
    auto [rhoLvals, rhoVvals] = get_funcvals(Theta, interval);
    // DCT to get the expansions in 1/m for each phase
    auto expL = get_expansion(rhoLvals, wmin, wmax);
    auto expV = get_expansion(rhoVvals, wmin, wmax);
    // Evaluate the 1D expansions
    return std::make_tuple(expL.y(1/m), expV.y(1/m));
}