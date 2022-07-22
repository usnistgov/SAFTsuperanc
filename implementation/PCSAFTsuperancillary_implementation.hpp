#pragma once

#include "PCSAFTsuperancillary_data.hpp"
#include "PCSAFTsuperancillary_crit.hpp"

using namespace ChebTools;
using namespace PCSAFTSuperAncillary;

const auto& get_interval(double w) {    
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
    return ChebyshevExpansion(
        std::move((V16 * vals.matrix()).eval().matrix()), wmin, wmax);
}

auto get_Ttilde_crit_min(double m){
    auto Ttilde_crit = cc_Ttilde(1/m);
    auto Ttilde_min = exp(-2.20078778)*pow(m, 0.37627892)*Ttilde_crit;
    return std::make_tuple(Ttilde_crit, Ttilde_min);
}

auto PCSAFTsuperanc_rhoLV(double Ttilde, double m){
    auto [Ttilde_crit, Ttilde_min] = get_Ttilde_crit_min(m);
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