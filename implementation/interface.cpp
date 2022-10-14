#include "PCSAFTsuperancillary_implementation.hpp"

#include "pcsaftsuperancversion.hpp"

#include "teqp/constants.hpp"

#if defined(PYBIND11) || defined(NANOBIND)

#if defined(PYBIND11)

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(PCSAFTsuperanc, m) {
    m.doc() = "SAFTsuperanc: Superancillary equations for the PC-SAFT EOS of Gross and Sadowski";
    m.attr("__version__") = PCSAFTSUPERANCVERSION;
    m.attr("N_A") = teqp::N_A;

    m.def("PCSAFTsuperanc_rhoLV", &PCSAFTsuperanc_rhoLV, py::arg("Ttilde"), py::arg("m"));
    m.def("get_Ttilde_crit_min", &get_Ttilde_crit_min, py::arg("m"));

}

#elif defined(NANOBIND)

#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/string.h>

namespace py = nanobind;

NB_MODULE(PCSAFTsuperanc, m) {
    m.attr("__doc__") = "SAFTsuperanc: Superancillary equations for the PC-SAFT EOS of Gross and Sadowski";
    m.attr("__version__") = PCSAFTSUPERANCVERSION;
    m.attr("N_A") = teqp::N_A;

    m.def("PCSAFTsuperanc_rhoLV", &PCSAFTsuperanc_rhoLV, py::arg("Ttilde"), py::arg("m"));
    m.def("get_Ttilde_crit_min", &get_Ttilde_crit_min, py::arg("m"));
}

#endif

#else 

int main(){
    double Ttilde = 0.462458913001258, m = 1.0;
    auto [rhotildeL, rhotildeV] = PCSAFTsuperanc_rhoLV(Ttilde, m);
    double rhoLexpected = 0.911529657714442, rhoVexpected = 0.0000169530410973542;
    int rr = 0; 
}

#endif 