#include <iostream>
#include "SuperAncillaryHelper.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>

#include "teqp/constants.hpp"
#include "teqp/models/pcsaft.hpp"
#include "teqp/algorithms/VLE.hpp"

#include "commons.hpp"

const std::string root = "../output";

TEST_CASE("Benchmark evaluation of critical ancillaries") {
	auto cch = CriticalCurveHelper(root);
	
	BENCHMARK("Call T") {
		return cch.Ttilde(1.0/5);
	};
	BENCHMARK("Call rho") {
		return cch.rhotilde(1.0/5);
	};
}

TEST_CASE("Benchmark matrix math") {
	Eigen::Matrix<double, 16, 16> Astatic16;
	Eigen::MatrixX<double> Adynamic16(16,16);
	Eigen::Matrix<double, 17, 17> Astatic17;
	Eigen::MatrixX<double> Adynamic17(17, 17);

	Eigen::Matrix<double, 16, 1> bstatic16;
	Eigen::MatrixX<double> bdynamic16(16, 1);
	Eigen::Matrix<double, 17, 1> bstatic17;
	Eigen::MatrixX<double> bdynamic17(17, 1);

	BENCHMARK("A*b (static) 16") {
		return Astatic16 * bstatic16;
	};
	BENCHMARK("A*b (static) 17") {
		return Astatic17 * bstatic17;
	};

	BENCHMARK("A*b (dynamic) 16") {
		return Adynamic16 * bdynamic16;
	};
	BENCHMARK("A*b (dynamic) 17") {
		return Adynamic17 * bdynamic17;
	};
}

TEST_CASE("Benchmark evaluation of super ancillaries", "[superanc]") {
	double w = 1 - 1.0/pow(2, 2);
	double mmin = 1.0, mmax = 1/(1*w + 1.0/64*(1 - w));
	auto anc = SuperAncillaryHelper<16>(root, mmin, mmax);
	auto cch = CriticalCurveHelper(root);

	double m = 5; 
	auto model5 = get_model(m);
	auto Ttildec = cch.Ttilde(1/m); 

	Eigen::VectorXd coef = Eigen::ArrayXd::LinSpaced(17, 0, 16);

	BENCHMARK("Get the values for nodes in w") {
		return anc.get_vals(0.8);
	}; 

	BENCHMARK("Get the expansions") {
		return anc.get_expansions(0.8);
	};

	BENCHMARK("Construct one ChebyshevExpansion") {
		return ChebTools::ChebyshevExpansion(coef, -5, 10);
	};
	
	BENCHMARK("Call single superancillary") {
		return anc(0.8, m);
	};

	BENCHMARK("Call superancillary and polish in double precision") {
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