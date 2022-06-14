#pragma once

#include "ChebTools/ChebTools.h"

// Data for a given interval in [wmin, wmax]
struct WInterval {
	const std::vector<double> wnodes;
	const std::vector<ChebTools::ChebyshevCollection> expsL, expsV; // In Theta, at each C-L node
	double wmin() const { return wnodes.front(); };
	double wmax() const { return wnodes.back(); };
};

struct CompleteWInterval {
	const std::vector<double> Wedges; 
	const std::vector<WInterval> intervals;
};

#include "domain.cpp"