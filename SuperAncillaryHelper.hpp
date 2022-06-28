# pragma once

#include <sstream>
#include <filesystem>
#include <fstream>

#include "nlohmann/json.hpp"

#include "ChebTools/ChebTools.h"

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
const Eigen::MatrixXd V8 = buildmat<8>();
const Eigen::MatrixXd V16 = buildmat<16>();

using namespace ChebTools;

/// Load a JSON file from a specified file
inline nlohmann::json load_JSON_file(const std::string& path) {
    if (!std::filesystem::is_regular_file(path)) {
        throw std::invalid_argument("Path to be loaded does not exist: " + path);
    }
    auto stream = std::ifstream(path);
    if (!stream) {
        throw std::invalid_argument("File stream cannot be opened from: " + path);
    }
    try {
        return nlohmann::json::parse(stream);
    }
    catch (...) {
        throw std::invalid_argument("File at " + path + " is not valid JSON");
    }
}

class CriticalCurveHelper {
private:
    auto get_coll(const std::string& root, const std::string& key) {
        std::string filepath = root + "/PCSAFT_crit_pts_expansions.json";
        auto j = load_JSON_file(filepath); 
        std::vector<ChebyshevExpansion> X;
        for (auto ex : j[key]) {
            auto get_array = [&ex](const std::string& k) {auto v = ex.at(k).get<std::vector<double>>();  return Eigen::Map<Eigen::ArrayXd>(&(v[0]), v.size()).eval(); };
            X.emplace_back(ChebyshevExpansion(get_array("coef"), ex.at("xmin"), ex.at("xmax")));
            double xmin = ex.at("xmin"), xmax = ex.at("xmax");
        }
        return ChebyshevCollection(X);
    }
public:
    ChebTools::ChebyshevCollection Ttilde, rhotilde;
    CriticalCurveHelper(const std::string& root) : Ttilde(get_coll(root, "Ttilde")), rhotilde(get_coll(root, "rhotilde")) {};
};

template<int Nm>
Eigen::Array<double, Nm + 1, 1> get_mnodes(double mmin, double mmax) {
    using ArrayN = Eigen::Array<double, Nm + 1, 1>;
    double ymin = 1 / mmax, ymax = 1 / mmin; // y = 1/m
    const ArrayN nodes_n11 = cos(ArrayN::LinSpaced(Nm + 1, 0, Nm)*EIGEN_PI/Nm); // [0,1,...,Nm]
    const ArrayN ynodes = (ymax - ymin)/2*nodes_n11 + (ymin + ymax) / 2;
    const ArrayN mnodes = 1/ynodes;
    return mnodes;
}

template<int Nm>
class SuperAncillaryHelper {

public:
    using ArrayN = Eigen::Array<double, Nm+1, 1>;
    const double mmin, mmax;
    const ArrayN mnodes;
    std::vector<ChebTools::ChebyshevCollection> expsL, expsV;
    Eigen::Matrix<double, Nm + 1, 2> F;
public:
    /// Unpack the expansions
    SuperAncillaryHelper(const std::string& root, double mmin, double mmax) : mmin(mmin), mmax(mmax), mnodes(get_mnodes<Nm>(mmin, mmax)){
        for (auto m : mnodes) {
            std::vector<ChebyshevExpansion> L, V;
            std::ostringstream ss; ss << std::scientific << std::setprecision(12) << root << "/PCSAFT_VLE_m" << m << "_expansions.json";
            std::string filepath = ss.str();
            if (!std::filesystem::exists(filepath)) {
                std::cout << mmin << " **** " << mmax << std::endl;
                std::cout << mnodes << std::endl;
                std::cout << "Missing file: " << filepath << std::endl;
            }
            for (auto ex : load_JSON_file(filepath)) {
                auto get_array = [&ex](const std::string &k) {auto v = ex.at(k).get<std::vector<double>>();  return Eigen::Map<Eigen::ArrayXd>(&(v[0]), v.size()).eval(); };
                Eigen::ArrayXd coefL = get_array("coefL"), coefV = get_array("coefV");
                L.emplace_back(ChebyshevExpansion(coefL, ex.at("xmin"), ex.at("xmax")));
                V.emplace_back(ChebyshevExpansion(coefV, ex.at("xmin"), ex.at("xmax")));
            }
            expsL.push_back(ChebyshevCollection(L));
            expsV.push_back(ChebyshevCollection(V));
        }
    }

    /// Call the function to get nodal values in 1/m
    auto get_vals(double Theta) const {
        ArrayN rhoLfvals, rhoVfvals;
        for (auto i = 0; i <= Nm; ++i) {
            rhoLfvals[i] = expsL[i](Theta);
            rhoVfvals[i] = expsV[i](Theta);
        }
        return std::make_tuple(rhoLfvals, rhoVfvals);
    }

    /// Retrieve a given ChebyshevCollection in Theta for a specified value of 1/m
    auto get_fitted_expansions(double m) const {
        for (auto i = 0; i < mnodes.size(); ++i) {
            if (mnodes[i] == m) {
                return std::make_tuple(expsL[i], expsV[i]);
            }
        }
        throw std::invalid_argument("Can't match m of " + std::to_string(m));
    }

    auto get_expansions(double Theta) {
        auto [rhoLfvals, rhoVfvals] = get_vals(Theta);
        double ymin = 1 / mmax, ymax = 1 / mmin;
        F.col(0) = rhoLfvals;
        F.col(1) = rhoVfvals;
        auto c = V16 * F;
        return std::make_tuple(
            ChebyshevExpansion(c.col(0).eval(), ymin, ymax),
            ChebyshevExpansion(c.col(1).eval(), ymin, ymax)
        );
    }

    /// Call the function to get densities from the superancillary
    auto operator()(double Theta, double m) {
        auto [expL, expV] = get_expansions(Theta);
        return std::make_tuple(expL.y(1/m), expV.y(1/m));
    }
};