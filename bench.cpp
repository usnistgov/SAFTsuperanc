#include "SuperAncillaryHelper.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>

#include "teqp/constants.hpp"
#include "teqp/models/pcsaft.hpp"
#include "teqp/algorithms/VLE.hpp"

#include "commons.hpp"

TEST_CASE("Benchmark evaluation of critical ancillaries") {
	auto cch = CriticalCurveHelper(".");
	
	BENCHMARK("Call T") {
		return cch.Ttilde(1.0/5);
	};
	BENCHMARK("Call rho") {
		return cch.rhotilde(1.0/5);
	};
}

TEST_CASE("Benchmark evaluation of super ancillaries") {
	auto anc = SuperAncillaryHelper<32>(".");
	auto cch = CriticalCurveHelper(".");

	double m = 5; 
	auto model5 = get_model(m);
	auto Ttildec = cch.Ttilde(1/m); 

	BENCHMARK("Get the values for nodes in 1/m") {
		return anc.getvals(0.8, m);
	}; 
	
	BENCHMARK("Call superancillary") {
		return anc(0.5, m);
	};

	BENCHMARK("Call superancillary and polish") {
		double T = 200, Ttilde = T/epsilon_over_k_K;
		std::valarray<double> c = { 0.37627892, -2.20078778 };
		auto Ttildemin = exp(c[1]) * pow(m, c[0])*Ttildec;
		double Theta = (Ttilde-Ttildemin)/(Ttildec - Ttildemin);
		auto [rhotildeL, rhotildeV] = anc(Theta, m);
		double rhoL = rhotildeL / (N_A * pow(sigma_m, 3));
		double rhoV = rhotildeV / (N_A * pow(sigma_m, 3));
		return teqp::pure_VLE_T(model5, T, rhoL, rhoV, 10);
	};
}