#include <sstream>

// Imports from boost
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
using namespace boost::multiprecision; 

#include "teqp/constants.hpp"
#include "teqp/derivs.hpp"
#include "teqp/models/pcsaft.hpp"
#include "teqp/algorithms/VLE.hpp"

#include "ChebTools/ChebTools.h"

// Global values, they cancel out; pick whatever you like
const double sigma_m = 1.0e-10;
const double epsilon_over_k_K = 100;
const double N_A = 6.022e23; // More digits are not important, they cancel

/// Build a PC-SAFT model for given value of m for a pure fluid. All other values are placeholders
inline auto get_model(const double m) {
    using namespace teqp::PCSAFT;
    std::vector<SAFTCoeffs> coeffs;
    SAFTCoeffs c;
    c.name = "PLACEHOLDER";
    c.m = m;
    c.sigma_Angstrom = sigma_m*1e10;
    c.epsilon_over_k = epsilon_over_k_K;
    c.BibTeXKey = "PLACEHOLDER";
    coeffs.push_back(c);
    return PCSAFTMixture(coeffs);
};

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
        return teqp::solve_pure_critical(get_model(m), Ttildec * epsilon_over_k_K, rhotildec / (N_A * pow(sigma_m, 3)));
    }
};

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 16)
{
    std::ostringstream out; out.precision(n);
    out << std::scientific << a_value;
    return out.str();
}

/// Trace away from the critical point down to a very low temperature
class VLETracer {
public:
    using my_float = double;
    using my_float_mp = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<200>>; 
    //using my_float_mp = double;
    std::vector<double> Ttilde, rhotildeL, rhotildeV;
    double Ttildec, rhotildec, m;

    VLETracer(my_float m, my_float Tc, my_float rhoc, double Tred_min) : m(m) {

        Ttildec = Tc / epsilon_over_k_K;
        rhotildec = rhoc * N_A * pow(sigma_m, 3);
        
        auto model = get_model(static_cast<double>(m));

        /// Store critical values are stored as tilde-scaled quantities
        Ttilde.push_back(static_cast<double>(Tc) / epsilon_over_k_K);
        rhotildeL.push_back(static_cast<double>(rhoc) * N_A * pow(sigma_m, 3));
        rhotildeV.push_back(static_cast<double>(rhoc) * N_A * pow(sigma_m, 3));

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
            const auto numeps = std::numeric_limits<my_float_mp>::epsilon();
            rhoL = rhovec[0], rhoV = rhovec[1];
            using std::isfinite;
            if (!isfinite(rhoL)) {
                break;
            }
            if (rhoV / rhoL < numeps) {
                break;
            }
            std::cout << T << " " << rhoL << " " << rhoV << " " << T/Tc << std::endl;
            // Values are stored as tilde-scaled quantities
            Ttilde.push_back(static_cast<double>(T) / epsilon_over_k_K);
            rhotildeL.push_back(static_cast<double>(rhoL) * N_A * pow(sigma_m, 3));
            rhotildeV.push_back(static_cast<double>(rhoV) * N_A * pow(sigma_m, 3));

            // Move down in temperature
            T -= 0.001*Tc;
            iter++;
        }
    }
    void to_json(const std::string path) {
        nlohmann::json outputs = {
            {"Ttilde", Ttilde},
            {"rhotildeL", rhotildeL},
            {"rhotildeV", rhotildeV},
            {"m", m},
            {"Ttildec", Ttildec},
            {"rhotildec", rhotildec},
        };
        std::ofstream file(path); file << outputs;
    }
    /// Linearly interpolate to get the tilde-scaled temperature and density
    auto interpolate_tilde(const double Ttilde_) {
        constexpr auto eps = std::numeric_limits<double>::epsilon();
        // Ttilde vector is stored in decreasing value, starting at critical and going down
        if (Ttilde_ > Ttilde[0]*(1+eps) || Ttilde_ < Ttilde[Ttilde.size() - 1]*(1-eps)) {
            throw std::invalid_argument("Ttilde is out of range");
        }
        else {
            for (auto i = 0; i < Ttilde.size() - 1; ++i) {
                if (Ttilde_ <= Ttilde[i]*(1+100*eps) && Ttilde_ >= Ttilde[i+1]*(1-100*eps)) {
                    return std::make_tuple(linterp(Ttilde, rhotildeL, i, Ttilde_), linterp(Ttilde, rhotildeV, i, Ttilde_));
                }
            }
            throw std::invalid_argument("How did I get here? Value of Ttilde is " + to_string_with_precision(Ttilde_) + " and limits are [" + to_string_with_precision(Ttilde[0]) + "," + to_string_with_precision(Ttilde[1]) + "]");
        }
    }

    auto build_expansions() {
        auto model = get_model(m);
        int Mnorm = 3;
        int max_refine_passes = 10; // As many as 2^max_refine_passes at end
        int Q = 0; // Vapor quality, also the index to rhovec to be returned
        auto f = [&](double Tred) {
            // Temperature to real units
            auto Ttilde = Tred * this->Ttildec;
            auto T = Ttilde * epsilon_over_k_K;
            
            // Interpolate to get approximate solution for tilde-scaled densities
            auto [rhotildeL, rhotildeV] = interpolate_tilde(Ttilde);
            auto rhoL = rhotildeL / (N_A * pow(sigma_m, 3));
            auto rhoV = rhotildeV / (N_A * pow(sigma_m, 3));

            // Solve phase equilibria problem in extended precision
            auto rhovec = teqp::pure_VLE_T<decltype(model), my_float_mp, teqp::ADBackends::multicomplex>(model, T, rhoL, rhoV, 10).cast<double>();
            std::cout << Tred << "::" << rhovec << std::endl;
            // Return the desired solution
            return rhovec[Q];
        };
        auto ces = ChebTools::ChebyshevExpansion::dyadic_splitting(12, f, Ttilde[0]/Ttildec, Ttilde.back()/Ttildec, Mnorm, 1e-12, max_refine_passes);
        for (auto ce : ces) {
            std::cout << "(" << ce.xmin() << "," << ce.xmax() << "): {" << ce.coef() << "}" << std::endl;
        }
    }
};

int main(){
    
    // Build interpolation function for critical point temperature and density as a function of 1/m
    double Ttilde1 = 1.2757487256161069; // nondimensional temperature
    double rhotilde1 = 0.2823886223463809; // nondimensional density
    Eigen::ArrayXd Ms = Eigen::ArrayXd::LinSpaced(100, 1, 100);
    auto ccrough = CriticalPointsInterpolator(Ms, Ttilde1*epsilon_over_k_K, rhotilde1/(N_A*pow(sigma_m, 3)));
    ccrough.to_json("PCSAFT_crit_pts_interpolation.json");

    for (double m : { 1,2,4,8,16,32,64 }) {
        // Solve for the exact critical point for this value of m
        auto [Tc, rhoc] = ccrough.get_Tcrhoc(m);
        // Trace to build a vector of values for saturation densities
        auto interp = VLETracer(m, Tc, rhoc, 0.001);
        interp.to_json("PCSAFT_VLE_m" + std::to_string(m) + ".json");
        try {
            interp.build_expansions();
        }
        catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
        }
    }
}