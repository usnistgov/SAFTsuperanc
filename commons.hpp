#pragma once

#include "teqp/models/pcsaft.hpp"
#include "teqp/algorithms/VLE.hpp"

// Global values, they cancel out; pick whatever you like
const double sigma_m = 1.0e-10;
const double epsilon_over_k_K = 100;
const double N_A = 6.022e23; // More digits are not important, they cancel

const auto dblepsilon = std::numeric_limits<double>::epsilon();

/// Build a PC-SAFT model for given value of m for a pure fluid. All other values are placeholders
inline auto get_model(const double m) {
	using namespace teqp::PCSAFT;
	std::vector<SAFTCoeffs> coeffs;
	SAFTCoeffs c;
	c.name = "PLACEHOLDER";
	c.m = m;
	c.sigma_Angstrom = sigma_m * 1e10;
	c.epsilon_over_k = epsilon_over_k_K;
	c.BibTeXKey = "PLACEHOLDER";
	coeffs.push_back(c);
	return PCSAFTMixture(coeffs);
};