/**
Test the obtained functions for their accuracy and dump results in JSON format for
further analysis in Python
*/

#include <iostream>
#include "SuperAncillaryHelper.hpp"
#include "nlohmann/json.hpp"

#include "commons.hpp"

void test_critical() {
	auto cch = CriticalCurveHelper(".");
	auto j = load_JSON_file("./PCSAFT_crit_pts_interpolation.json");
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

template<int Nm>
void test_VLE(const Eigen::ArrayXd& mvec, const std::string& path) {
	auto cch = CriticalCurveHelper(".");
	auto anc = SuperAncillaryHelper<Nm>(".");
	std::vector<nlohmann::json> jout;
	int Ntheta = 100;
	double dTheta = 1.0 / (Ntheta - 1.0);
	for (double m : mvec) {
		std::cout << "m: " << m << std::endl;
		auto model = get_model(m);
		double epsilon_over_k_K = model.get_epsilon_over_k_K()[0];
		double sigma_m = model.get_sigma_Angstrom()[0]/1e10;
		double sigma_m_3 = pow(sigma_m, 3);
		for (auto Theta = 0.0; Theta <= 1; Theta += dTheta) {
			// Unpack Theta to practical temperature variables
			auto Ttildec = cch.Ttilde(1/m);
			auto Ttilde_min = exp(-2.20078778) * pow(m, 0.37627892) * Ttildec;
			auto Ttilde = Ttilde_min + Theta * (Ttildec - Ttilde_min);
			auto T = Ttilde*epsilon_over_k_K;

			// Evaluate the ancillary functions
			auto [rhotildeL, rhotildeV] = anc(Theta, m);
			double rhoL = rhotildeL / (N_A * sigma_m_3);
			double rhoV = rhotildeV / (N_A * sigma_m_3);

			// Polish the solution, in double and then extended precision
			using my_float_mp = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<200>>;
			Eigen::ArrayXd rhovec_dbl = teqp::pure_VLE_T(model, T, rhoL, rhoV, 10);
			Eigen::ArrayXd rhovec_mp = teqp::pure_VLE_T<decltype(model), my_float_mp, teqp::ADBackends::multicomplex>(model, T, rhovec_dbl[0], rhovec_dbl[1], 10).cast<double>();
			double rhoLerr = rhovec_mp[0]*N_A*sigma_m_3/rhotildeL - 1;
			std::cout << rhoLerr << std::endl;
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
	}
	std::ofstream file(path); file << jout;
}

int main() {
	//test_critical();
	double mmin = 1, mmax = 64;
	//test_VLE<32>(get_mnodes<32>(1, 64), "PCSAFT_VLE_check_fitted.json");
	
	// For an expansion of degree N, the worst-case error should be at the Chebyshev-Lobatto nodes of degree 2N
	// These are the odd nodes of the expansion of degree 2N
	auto worst_case_nodes = [](const Eigen::ArrayXd& nodes) {return nodes(Eigen::seq(1, Eigen::last, 2)).eval(); };
	test_VLE<64>(worst_case_nodes(get_mnodes<128>(1, 64)), "PCSAFT_VLE_check_worstcase_Nm64.json"); 
	test_VLE<32>(worst_case_nodes(get_mnodes<64>(1, 64)),  "PCSAFT_VLE_check_worstcase_Nm32.json");
	test_VLE<16>(worst_case_nodes(get_mnodes<32>(1, 64)),  "PCSAFT_VLE_check_worstcase_Nm16.json");
}