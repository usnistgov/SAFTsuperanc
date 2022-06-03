#include <sstream>
#include <filesystem>

// Imports from boost
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost::multiprecision; 
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "teqp/constants.hpp"
#include "teqp/derivs.hpp"

#include "ChebTools/ChebTools.h"

#include "commons.hpp"



template<typename T>
inline auto linterp (const T& x, const T& y, std::size_t i, const double val) {
    return (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (val - x[i]) + y[i];
};

/// A class for building an interpolator of critical points as a function of 1/m
/// 
/// Arrays of temperatures and densities are in tilde-reduced form
template<typename Ms>
class CriticalPointsInterpolator {
private:
    const Ms ms;
    std::vector<double> Ttilde, rhotilde, one_over_m;
public:
    CriticalPointsInterpolator(const Ms & ms, double T0, double rho0) : ms(ms) {
        double T = T0, rho = rho0;
        for (double m : ms) {
            auto model = get_model(m);
            std::tie(T, rho) = teqp::solve_pure_critical(model, T, rho);
            if (!std::isfinite(T)){
                throw std::invalid_argument("Critical point solving failed for m of " + std::to_string(m));
            }
            Ttilde.push_back(T/epsilon_over_k_K);
            rhotilde.push_back(rho*N_A*pow(sigma_m, 3));
            one_over_m.push_back(1/m);
            std::cout << "[crit] " <<  m << "::" << T << "," << rho << std::endl;
        }
    };
    void to_json(const std::string path) {
        nlohmann::json outputs = {
            {"Ttilde", Ttilde},
            {"rhotilde", rhotilde},
            {"1/m", one_over_m},
        };
        std::ofstream file(path);
        file << outputs;
    }
    /// Linearly interpolate to get the tilde-scaled temperature and density
    auto interpolate_tilde(const double m) {
        if (m < ms[0] || m > ms[ms.size()-1]) {
            throw std::invalid_argument("m is out of range");
        }
        else {
            for (auto i = 0; i < ms.size()-1; ++i) {
                if (m >= ms[i] && m <= ms[i + 1]) {
                    return std::make_tuple(linterp(one_over_m, Ttilde, i, 1/m), linterp(one_over_m, rhotilde, i, 1/m));
                }
            }
            throw std::invalid_argument("How did I get here for m of " + std::to_string(m) + "for critical point interpolation");
        }
    }
    auto get_Tcrhoc(double m) {
        auto [Ttildec, rhotildec] = interpolate_tilde(m);
        return teqp::solve_pure_critical(get_model(m), Ttildec*epsilon_over_k_K, rhotildec / (N_A * pow(sigma_m, 3)));
    }

    auto build_expansions(double mmin, double mmax) {
        double ymin = 1 / mmax, ymax = 1 / mmin;
        auto f_ = [this](double y, int i) {
            double m = 1/y;
            auto [T, rho] = get_Tcrhoc(m);
            if (!std::isfinite(T)) {
                throw std::invalid_argument("Critical point solving failed for m of " + std::to_string(m));
            }
            if (i == 0) {
                return T/epsilon_over_k_K;
            }
            else {
                return rho*N_A*pow(sigma_m, 3);
            }
        };
        auto f_T = [f_](double y) { return f_(y, 0); };
        auto f_rho = [f_](double y) { return f_(y, 1); };

        int Mnorm = 3;
        double tol = 1e-12;
        int max_split = 12;
        return std::make_tuple(
            ChebTools::ChebyshevExpansion::dyadic_splitting<std::vector<ChebTools::ChebyshevExpansion>>(16, f_T, ymin, ymax, Mnorm, tol, max_split),
            ChebTools::ChebyshevExpansion::dyadic_splitting<std::vector<ChebTools::ChebyshevExpansion>>(16, f_rho, ymin, ymax, Mnorm, tol, max_split)
        );
    }
    auto dump_expansions(double mmin, double mmax, const std::string& path) {
        auto [exsT, exsrho] = build_expansions(mmin, mmax);
        std::vector<std::string> keys = { "Ttilde", "rhotilde" };
        std::map<std::string, nlohmann::json> outputs;
        for (auto i = 0; i < keys.size(); ++i) {
            auto expansions = (i == 0) ? exsT : exsrho;
            std::vector<nlohmann::json> jexpansions;
            for (auto ex : expansions) {
                jexpansions.push_back({
                    {"xmin", ex.xmin()},
                    {"xmax", ex.xmax()},
                    {"coef", ex.coef()},
                    });
            }
            outputs[keys[i]] = jexpansions;
        }
        std::ofstream(path) << outputs;
    }
};

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 16)
{
    std::ostringstream out; out.precision(n);
    out << std::scientific << a_value;
    return out.str();
}

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

/// Trace away from the critical point down to a very low temperature
class VLETracer {
public:
    using my_float = double;
    using my_float_mp = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<200>>; 
    std::vector<double> Ttilde, rhotildeL, rhotildeV;
    double Ttildec, rhotildec, m;

    VLETracer(my_float m, my_float Tc, my_float rhoc, double Tred_min) : m(m) {

        Ttildec = Tc / epsilon_over_k_K;
        rhotildec = rhoc * N_A * pow(sigma_m, 3);
        
        auto model = get_model(static_cast<double>(m));

        /// Store critical values are stored as tilde-scaled quantities
        Ttilde.push_back(static_cast<double>(Tc) / epsilon_over_k_K);
        rhotildeL.push_back(static_cast<double>(rhoc)*N_A*pow(sigma_m, 3));
        rhotildeV.push_back(static_cast<double>(rhoc)*N_A*pow(sigma_m, 3));

        // Move just away from the critical point
        // with the critical point extrapolation
        my_float T = Tc*(1 - 1e-3);
        auto rhos = teqp::extrapolate_from_critical(model, Tc, rhoc, T);

        // Working variables are in molar units
        my_float_mp rhoL = rhos[0], rhoV = rhos[1];
        int iter = 0;
        while (T > Tred_min*Tc) {
            
            // Do the VLE calculation
            //auto rhovec = teqp::pure_VLE_T(model, T, rhoL, rhoV, 10);
            auto rhovec = teqp::pure_VLE_T<decltype(model), my_float_mp, teqp::ADBackends::multicomplex>(model, T, rhoL, rhoV, 10);

            teqp::IsothermPureVLEResiduals<decltype(model), my_float_mp, teqp::ADBackends::multicomplex> res(model, T);
            auto residuals = res.call(rhovec);
            if (residuals.abs().maxCoeff() > 1e-16) {
                std::cout << "Large residual detected : " << residuals.maxCoeff() << std::endl;
            }

            const auto numeps = std::numeric_limits<my_float_mp>::epsilon();
            rhoL = rhovec[0], rhoV = rhovec[1];
            using std::isfinite;
            if (!isfinite(rhoL)) {
                std::cout << "not finite" << std::endl;
                break;
            }
            if (rhoV / rhoL < numeps) {
                break;
            }
            auto max_rho = model.max_rhoN(T, (Eigen::ArrayXd(1) << 1.0).finished())/N_A;
            std::cout << T << " " << rhoL << " " << rhoV << " " << T/Tc <<  " " << rhoL/max_rho << std::endl;
            // Values are stored as tilde-scaled quantities
            Ttilde.push_back(static_cast<double>(T) / epsilon_over_k_K);
            rhotildeL.push_back(static_cast<double>(rhoL)*N_A*pow(sigma_m, 3));
            rhotildeV.push_back(static_cast<double>(rhoV)*N_A*pow(sigma_m, 3));

            // Move down in temperature
            T -= 0.00025*Tc;
            iter++;
        }
    }
    VLETracer(const std::string &path) {
        auto j = load_JSON_file(path);
        m = j.at("m");
        Ttildec = j.at("Ttildec");
        rhotildec = j.at("rhotildec"); 
        Ttilde = j.at("Ttilde").get<std::vector<double>>();
        rhotildeL = j.at("rhotildeL").get<std::vector<double>>();
        rhotildeV = j.at("rhotildeV").get<std::vector<double>>();
    }
    void to_json(const std::string path) {
        nlohmann::json outputs = {
            {"m", m},
            {"Ttildec", Ttildec},
            {"rhotildec", rhotildec},
            {"Ttilde", Ttilde},
            {"rhotildeL", rhotildeL},
            {"rhotildeV", rhotildeV},
        };
        std::ofstream file(path); file << outputs;
    }

    /// Linearly interpolate to get the tilde-scaled temperature and density
    auto interpolate_tilde(const double Ttilde_) {
        constexpr auto eps = std::numeric_limits<double>::epsilon();
        // Ttilde vector is stored in decreasing value, starting at critical and going down
        if (Ttilde_ > Ttilde[0]*(1+100*eps) || Ttilde_ < Ttilde.back()*(1-100*eps)) {
            throw std::invalid_argument("Ttilde is out of range for value of " + to_string_with_precision(Ttilde_) + " and limits are [" + to_string_with_precision(Ttilde[0]) + "," + to_string_with_precision(Ttilde[Ttilde.size() - 1]) + "]");
        }
        else {
            for (auto i = 0; i < Ttilde.size() - 1; ++i) {
                if (Ttilde_ <= Ttilde[i]*(1+100*eps) && Ttilde_ >= Ttilde[i+1]*(1-100*eps)) {
                    return std::make_tuple(linterp(Ttilde, rhotildeL, i, Ttilde_), linterp(Ttilde, rhotildeV, i, Ttilde_));
                }
            }
            throw std::invalid_argument("How did I get here? Value of Ttilde is " + to_string_with_precision(Ttilde_) + " and limits are [" + to_string_with_precision(Ttilde[0]) + "," + to_string_with_precision(Ttilde[Ttilde.size() - 1]) + "]");
        }
    }

    auto test_interpolation(const std::size_t N) {
        auto r01 = (Eigen::ArrayXd::Random(N)/2+0.5).eval();
        auto model = get_model(m);
        auto r01min = r01.minCoeff(), r01max = r01.maxCoeff();
        auto Ttildevals = r01*Ttildec + (1.0-r01)*Ttilde.back();
        auto Ttmin = Ttildevals.minCoeff(), Ttmax = Ttildevals.maxCoeff();
        for (double Ttilde : Ttildevals) {
            // Do the VLE calculation
            auto [rhoLtilde, rhoVtilde] = interpolate_tilde(Ttilde);
            std::cout << "." << std::endl;
            auto T = Ttilde*epsilon_over_k_K;
            auto rhovec = teqp::pure_VLE_T<decltype(model), my_float_mp, teqp::ADBackends::multicomplex>(model, T, rhoLtilde/(N_A * pow(sigma_m, 3)), rhoVtilde/(N_A * pow(sigma_m, 3)), 10);
            using std::isfinite;
            if (!isfinite(rhovec[0])) {
                std::cout << Ttilde/Ttildec << std::endl;
            }
        }
    }

    auto test_extrapolation_from_critical(const double m, double Tc, double rhoc) {
        auto model = get_model(m);
        for (double dT = 0.01*Tc; dT > 1e-10*Tc; dT /= 2){
            // Working variables are in molar units, as is extrapolation because it uses teqp
            auto rhoe = teqp::extrapolate_from_critical(model, Tc, rhoc, Tc-dT);
            // Do the VLE calculation in extended precision
            auto rhovec = teqp::pure_VLE_T<decltype(model), my_float_mp, teqp::ADBackends::multicomplex>(model, Tc-dT, rhoe[0], rhoe[1], 10);
            std::cout << dT << " ;; " << rhovec << std::endl;
        }
    }

    /// Get value of Ttilde for which density ratio is approximately 10^20
    double get_Ttildemin(double m) {
        std::valarray<double> c = { 0.37627892, -2.20078778 };
        auto Tred = exp(c[1])*pow(m, c[0]);
        return Tred * Ttildec;
    }

    /// Build all the expansions
    auto build_expansions(const std::string &path) {
        if (std::filesystem::exists(path)) {
            return;
        }
        auto model = get_model(m);
        int Mnorm = 3;
        // Convenience function to get the M-element norm
        auto get_err = [Mnorm](const auto& ce) { return ce.coef().tail(Mnorm).norm() / ce.coef().head(Mnorm).norm(); };
        int max_refine_passes = 10; // As many as 2^max_refine_passes at end
        
        // The data type to contain the two expansions
        struct ChebPair { ChebTools::ChebyshevExpansion rhoL, rhoV; };
        std::vector<ChebPair> expansions;
        auto Ttildemin = get_Ttildemin(m);

        /// A function to return liquid and vapor densities given temperature
        /// in normalized form (Theta = (T-Tmin)/(Tc-Tmin))
        /// Returned densities are in tilde-reduced form as an Eigen::Array
        auto do_VLE = [&expansions, &model, &Ttildemin, this](double Theta) {

            // From reduced scaled temperature to real units
            auto Ttilde = Theta * (Ttildec - Ttildemin) + Ttildemin;
            auto Tred = Ttilde / Ttildec;
            auto T = Ttilde * epsilon_over_k_K;

            double rhotildeL, rhotildeV;
            bool interpolation_used = false;
            if (std::abs(Tred - 1) < dblepsilon * 10) {
                // At the critical point (to numerical precision) 
                // There can be a small loss in precision caused by intermediate calculations such that the limit for Theta may be no longer precisely 1.0
                return (Eigen::Array<double, 2, 1>() << rhotildec,rhotildec).finished();
            }
            else if (Ttilde > this->Ttilde[1]) {
                auto rhos = teqp::extrapolate_from_critical(model, Ttildec*epsilon_over_k_K, rhotildec/(N_A*pow(sigma_m, 3)), T);
                rhotildeL = rhos[0]*(N_A*pow(sigma_m, 3));
                rhotildeV = rhos[1]*(N_A*pow(sigma_m, 3));
            }
            else{
                // Interpolate from the rough values to get approximate solution for tilde-scaled densities
                std::tie(rhotildeL, rhotildeV) = interpolate_tilde(Ttilde);
                interpolation_used = true;
            }
            auto rhoL = rhotildeL/(N_A*pow(sigma_m, 3));
            auto rhoV = rhotildeV/(N_A*pow(sigma_m, 3));

            // Solve phase equilibria problem in extended precision
            auto rhovec = teqp::pure_VLE_T<decltype(model), my_float_mp, teqp::ADBackends::multicomplex>(model, T, rhoL, rhoV, 10).cast<double>();
            ////std::cout << Tred << "::" << rhovec << std::endl;
            if (!std::isfinite(rhovec[0])) {
                throw std::invalid_argument("Invalid VLE solution for starting values of (" + to_string_with_precision(rhoL) 
                    + "," + to_string_with_precision(rhoV) + ") mol/m^3 for Tred of " + to_string_with_precision(Tred) + ((interpolation_used) ? " [interpolation]" : " [critical]"));
            }
            // Return the desired solutions in tilde-reduced form
            return (rhovec* (N_A * pow(sigma_m, 3))).eval();
        };

        auto make_expansions = [&](const std::size_t N, double Tmin, double Tmax) -> ChebPair {
            // Return liquid and vapor expansions
            // 
            // Node values in [-1, 1]
            auto nodesn11 = ChebTools::get_CLnodes(N).array();
            auto nodes = ((Tmax - Tmin) * nodesn11 + (Tmax + Tmin)) / 2;
            // Densities of both phases at node values
            Eigen::ArrayXd rhoLvals(nodes.size()), rhoVvals(nodes.size());
            for (auto i = 0; i < nodes.size(); ++i) {
                auto rhos = do_VLE(nodes[i]);
                rhoLvals[i] = rhos[0];
                rhoVvals[i] = rhos[1];
            }
            //std::cout << rhoLvals << std::endl;
            //std::cout << rhoVvals << std::endl;
            // Invert values with DCT to get coefficients and return the expansions
            return ChebPair{ 
                ChebTools::ChebyshevExpansion::factoryf(N, rhoLvals, Tmin, Tmax), 
                ChebTools::ChebyshevExpansion::factoryf(N, rhoVvals, Tmin, Tmax) };
        };
        
        int N = 12;
        double tol = 1e-12;
        double Theta_min = 0.0, Theta_max = 1.0;
        expansions.emplace_back(make_expansions(N, Theta_min, Theta_max));

        // Now enter into refinement passes
        for (int refine_pass = 0; refine_pass < max_refine_passes; ++refine_pass) {
            std::cout << "pass " << refine_pass << std::endl;
            bool all_converged = true;
            // Start at the right and move left because insertions will make the length increase
            for (int iexpansion = static_cast<int>(expansions.size()) - 1; iexpansion >= 0; --iexpansion) {
                auto& expan = expansions[iexpansion];
                auto errL = get_err(expan.rhoL);
                auto errV = get_err(expan.rhoV);
                if (errL > tol || errV > tol) {
                    // Splitting is required, do a dyadic split
                    auto xmid = (expan.rhoL.xmin() + expan.rhoL.xmax()) / 2;
                    auto newlefts = make_expansions(N, expan.rhoL.xmin(), xmid);
                    auto newrights = make_expansions(N, xmid, expan.rhoL.xmax());

                    // Function to check if any coefficients are invalid (evidence of a bad function value)
                    auto all_coeffs_ok = [](const auto& v) {
                        for (auto i = 0; i < v.size(); ++i) {
                            if (!std::isfinite(v[i])) { return false; }
                        }
                        return true;
                    };
                    // Check if any coefficients are invalid, stop if so
                    if (   !all_coeffs_ok(newlefts.rhoL.coef()) || !all_coeffs_ok(newrights.rhoL.coef()) 
                        || !all_coeffs_ok(newlefts.rhoV.coef()) || !all_coeffs_ok(newrights.rhoV.coef())) {
                        throw std::invalid_argument("At least one coefficient is non-finite");
                    }
                    std::swap(expan, newlefts);
                    expansions.insert(expansions.begin() + iexpansion + 1, newrights);
                    all_converged = false;
                }
            }
            if (all_converged) { break; }
        }        
     
        std::vector<nlohmann::json> jexpansions;
        for (auto ex : expansions) {
            std::cout << "(" << ex.rhoL.xmin() << "," << ex.rhoV.xmax() << "): {" << ex.rhoL.coef() << "}" << std::endl;
            jexpansions.push_back({
                {"xmin", ex.rhoL.xmin()},
                {"xmax", ex.rhoL.xmax()},
                {"coefL", ex.rhoL.coef()},
                {"coefV", ex.rhoV.coef()}
            });
        }
        std::ofstream file(path); file << jexpansions;
    }
};

int main(){
    
    // Build interpolation function for critical point temperature and density as a function of 1/m
    double Ttilde1 = 1.2757487256161069; // nondimensional temperature
    double rhotilde1 = 0.2823886223463809; // nondimensional density
    Eigen::ArrayXd Ms = Eigen::ArrayXd::LinSpaced(1000, 1, 1000);
    auto ccrough = CriticalPointsInterpolator(Ms, Ttilde1*epsilon_over_k_K, rhotilde1/(N_A*pow(sigma_m, 3)));
    ccrough.to_json("PCSAFT_crit_pts_interpolation.json");
    double mmincrit = 1.0, mmaxcrit = 900.0;
    ccrough.dump_expansions(mmincrit, mmaxcrit, "PCSAFT_crit_pts_expansions.json");

    // Chebyshev nodes in y=1/m for m from 1 to mmax
    auto N = 8;
    double mmin = 1.0, mmax = 32.0;
    double ymin = 1.0 / mmax, ymax = 1.0/mmin;
    auto nodesn11 = ChebTools::get_CLnodes(N).array();
    auto ynodes = ((ymax - ymin) * nodesn11 + (ymax + ymin)) / 2;
    auto mnodes = 1 / ynodes;

    boost::asio::thread_pool pool(4); // 4 threads
    
    for (double m : mnodes) {
        auto one_m = [m, &ccrough] {
            std::cout << "************** " << m << " *****************" << std::endl;

            // Solve for the exact critical point for this value of m
            auto [Tc, rhoc] = ccrough.get_Tcrhoc(m);
            // Trace to build a vector of values for saturation densities
            std::string path = "PCSAFT_VLE_m" + std::to_string(m) + ".json";
            auto interp = (std::filesystem::exists(path)) ? VLETracer(path) : VLETracer(m, Tc, rhoc, 0.001);
            interp.to_json(path);
            //interp.test_interpolation(10000);
            //interp.test_extrapolation_from_critical(m, Tc, rhoc);

            try {
                std::string path_expansions = "PCSAFT_VLE_m" + std::to_string(m) + "_expansions.json";
                interp.build_expansions(path_expansions);
            }
            catch (const std::exception& e) {
                std::cout << e.what() << std::endl;
            }
        };
        boost::asio::post(pool, one_m);
    }
    pool.join();
}