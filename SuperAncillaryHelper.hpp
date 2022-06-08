# pragma once

#include <sstream>
#include <filesystem>
#include <fstream>

#include "nlohmann/json.hpp"

#include "ChebTools/ChebTools.h"

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
    auto getvals(double Theta) {
        ArrayN rhoLfvals, rhoVfvals;
        for (auto i = 0; i <= Nm; ++i) {
            rhoLfvals[i] = expsL[i](Theta);
            rhoVfvals[i] = expsV[i](Theta);
        }
        return std::make_tuple(rhoLfvals, rhoVfvals);
    }

    auto get_expansions(double Theta) {
        auto [rhoLfvals, rhoVfvals] = getvals(Theta);
        double ymin = 1 / mmax, ymax = 1 / mmin;
        return std::make_tuple(
            ChebyshevExpansion::factoryf(Nm, rhoLfvals, ymin, ymax),
            ChebyshevExpansion::factoryf(Nm, rhoVfvals, ymin, ymax)
        );
    }

    /// Call the function to get densities from the superancillary
    auto operator()(double Theta, double m) {
        auto [rhoLfvals, rhoVfvals] = getvals(Theta);
        auto [expL, expV] = get_expansions(Theta);
        return std::make_tuple(expL.y(1/m), expV.y(1/m));
    }
};