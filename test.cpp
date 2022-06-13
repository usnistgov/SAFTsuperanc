/**
Test the obtained functions for their accuracy and dump results in JSON format for
further analysis in Python
*/

#include "teqp/types.hpp"

#include <iostream>
#include "SuperAncillaryHelper.hpp"
#include "nlohmann/json.hpp"

#include "commons.hpp"

const std::string root = "../output";

void test_critical() {
	auto cch = CriticalCurveHelper(root);
	auto j = load_JSON_file(root + "/PCSAFT_crit_pts_interpolation.json");
	std::vector<nlohmann::json> jout;
	auto N = j["Ttilde"].size();
	for (auto i = 0; i < N; ++i) {
		double mrecip = j["1/m"][i], m = 1/mrecip;
		if (m > 1/cch.Ttilde.get_exps()[0].xmin()) {
			continue;
		}
		double Ttilde_fit = cch.Ttilde(mrecip);
		double rhotilde_fit = cch.rhotilde(mrecip);
		double Ttilde_tab = j.at("Ttilde")[i];
		double rhotilde_tab = j.at("rhotilde")[i];
		jout.push_back({
			{"1/m", mrecip},
			{"Ttilde_fit", Ttilde_fit},
			{"Ttilde_tab", Ttilde_tab},
			{"rhotilde_fit", rhotilde_fit},
			{"rhotilde_tab", rhotilde_tab}
		});
	}
	std::ofstream file("PCSAFT_crit_pts_check.json"); file << jout;
}

enum class node_options { normal, worst };

template<int Nm, node_options opt>
void test_VLE(double mmin, double mmax, int m_split, const std::string& path) {
	auto cch = CriticalCurveHelper(root);
	auto Wedges = load_JSON_file(root + "/Wedges_pass" + std::to_string(m_split) + ".json").at("Wedges").get<std::vector<double>>();
	for (int i = static_cast<int>(Wedges.size())-2; i >= 0; --i) { 
		double ymin = Wedges[i], ymax = Wedges[i + 1];
		double mmin_ = 1 / ymax, mmax_ = 1 / ymin;
		std::cout << mmin_ << "," << mmax_ << std::endl;

		try{
		auto anc = SuperAncillaryHelper<Nm>(root, 1/ymax, 1/ymin);
		std::vector<nlohmann::json> jout;
		int Ntheta = 100;
		int j = 0;
		constexpr int Nm_domain = [&] {return (opt == node_options::normal) ? Nm : 2 * Nm; }();
		for (double m : get_mnodes<Nm_domain>(1/ymax, 1/ymin)) {
			if (opt == node_options::worst && j % 2 == 0) {
				// Only keep the odd nodes if doing worst case error calcs
				j++; continue;
			}
			
			//continue;

			auto model = get_model(m);
			double epsilon_over_k_K = model.get_epsilon_over_k_K()[0];
			double sigma_m = model.get_sigma_Angstrom()[0]/1e10;
			double sigma_m_3 = pow(sigma_m, 3);

			Eigen::ArrayXd Thetavec = Eigen::ArrayXd::LinSpaced(Ntheta + 1, 0, 1);
			if (opt == node_options::normal) {

				/// All the nodes of the refined expansion for this value of w
				auto [expL, expV] = anc.get_fitted_expansions(m);
				std::vector<double> Thetas;
				for (auto ki = 0; ki < expL.get_exps().size(); ++ki) { // Iterate over intervals in Theta
					double xmin = expL.get_exps()[ki].xmin(), xmax = expL.get_exps()[ki].xmax();
					int N = static_cast<int>(expL.get_exps()[ki].coef().size()) - 1;
					// Double-sample in the Theta direction for this expansion
					auto Thetavals = (xmax - xmin) / 2 * ChebTools::get_CLnodes(2 * N).array() + (xmax + xmin) / 2;
					// Odd nodes (worst-case error)
					for (auto Thetavecodd : Thetavals(Eigen::seq(1, 2 * N + 1, 2))) {
						Thetas.push_back(Thetavecodd);
					};
				}
				Thetavec = Eigen::Map<Eigen::ArrayXd>(&(Thetas[0]), Thetas.size());
				
			}
			std::cout << "m: " << m << "; len(Thetavec): " << Thetavec.size() << std::endl;

			for (auto Theta : Thetavec) {
				// Unpack Theta to practical temperature variables
				auto Ttildec = cch.Ttilde(1/m);
				auto Ttilde_min = exp(-2.20078778)*pow(m, 0.37627892)*Ttildec;
				auto Ttilde = Ttilde_min + Theta*(Ttildec - Ttilde_min);
				auto T = Ttilde*epsilon_over_k_K;

				// Evaluate the ancillary functions
				auto [rhotildeL, rhotildeV] = anc(Theta, m);
				double rhoL = rhotildeL/(N_A*sigma_m_3);
				double rhoV = rhotildeV/(N_A*sigma_m_3);

				// Polish the solution, in double and then extended precision
				using my_float_mp = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<50>>;
				Eigen::ArrayXd rhovec_dbl = teqp::pure_VLE_T(model, T, rhoL, rhoV, 10);
				//Eigen::ArrayXd rhovec_mp = -100 * rhovec_dbl;
				Eigen::ArrayXd rhovec_mp = teqp::pure_VLE_T<decltype(model), my_float_mp, teqp::ADBackends::multicomplex>(model, T, rhoL, rhoV, 10).cast<double>();
				double rhoLerr = rhovec_dbl[0]*N_A*sigma_m_3/rhotildeL - 1;

				if (std::abs(Theta - 1) < dblepsilon) {
					auto rhoc = cch.rhotilde(1/m)/(N_A*sigma_m_3);
					rhovec_dbl << rhoc, rhoc;
					rhovec_mp << rhoc, rhoc;
				}
				//std::cout << rhoLerr << std::endl;
				jout.push_back({
					{"1/m", 1/m},
					{"rhotildeL_anc", rhotildeL},
					{"rhotildeL_VLE", rhovec_dbl[0]*N_A*sigma_m_3},
					{"rhotildeL_VLEmp", rhovec_mp[0]*N_A*sigma_m_3},
					{"rhotildeV_anc", rhotildeV},
					{"rhotildeV_VLE", rhovec_dbl[1]*N_A*sigma_m_3},
					{"rhotildeV_VLEmp", rhovec_mp[1]*N_A*sigma_m_3},
					{"Theta", Theta},
					{"Ttilde", Ttilde},
					{"Ttilde_crit", Ttildec},
					{"Ttilde_min", Ttilde_min}
				});
				//std::cout << jout.back() << std::endl;
			}
			j++;
		}
		std::ofstream file(std::to_string(i) + path); file << jout;
		}catch(...){}
	}
}

int main() {
	test_critical();
    
	double mmin = 1.0, mmax = 64.0;
	int m_split = 8;

	// The nodes where the fitting was done, should be good to close to numerical precision
	test_VLE<16, node_options::normal>(mmin, mmax, m_split, "PCSAFT_VLE_check_fitted_Nm16.json");

	// For an expansion of degree N, the worst-case error should be at or close to the odd Chebyshev-Lobatto nodes of degree 2N
	// (the even ones should be good to numerical precision)
	//test_VLE<16, node_options::worst>(mmin, mmax, m_split, "PCSAFT_VLE_check_worstcase_Nm16.json");
}