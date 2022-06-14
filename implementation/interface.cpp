#include "PCSAFTsuperancillary_implementation.hpp"

#if defined(PYBIND11)

#include "pcsaftsuperancversion.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

void init_superanc(py::module& m) {
    m.def("PCSAFTsuperanc_rhoLV", &PCSAFTsuperanc_rhoLV, py::arg("Tilde"), py::arg("m"));
}

PYBIND11_MODULE(PCSAFTsuperanc, m) {
    m.doc() = "SAFTsuperanc: Superancillary equations for the PC-SAFT EOS of Gross and Sadowski";
    m.attr("__version__") = PCSAFTSUPERANCVERSION;
    init_superanc(m);
}

#else 

int main(){
    double Ttilde = 0.462458913001258, m = 1.0;
    auto [rhotildeL, rhotildeV] = PCSAFTsuperanc_rhoLV(Ttilde, m);
    double rhoLexpected = 0.911529657714442, rhoVexpected = 0.0000169530410973542;
    int rr = 0; 
}

#endif 