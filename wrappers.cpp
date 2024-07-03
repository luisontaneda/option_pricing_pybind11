#include <pybind11/pybind11.h>

#include "pricing_algorithms/COS_method/cos_pricing.h"
#include "pricing_algorithms/montecarlo_method/stoch_proccess.h"

namespace py = pybind11;

void init_montecarlo_sub(py::module &);
void init_vanilla_mont_sub(py::module &);
void init_european_mont_sub(py::module &);


void init_asian_cos_sub(py::module &m);
void init_vanilla_cos_sub(py::module &m);
void init_exotic_cos_sub(py::module &m);
void init_european_cos_sub(py::module &m);
void init_cos_method_sub(py::module &m);


PYBIND11_MODULE(option_pricing, m) {
    m.doc() = "option_pricing plugin"; // optional module docstring

    py::module montecarlo_sub = m.def_submodule("montecarlo", "A submodule for Monte Carlo methods");
    py::module cos_method_sub = m.def_submodule("cos_method", "A submodule for COS method");
    
    // Initialize the submodules
    init_montecarlo_sub(montecarlo_sub);
    init_cos_method_sub(cos_method_sub);
}


// submodules for Montecarlo method


void init_asian_mont_sub(py::module &m) {
    m.def("gbm", &asian_gbm, "Pricing using Stochastic process: GBM");
}


void init_vanilla_mont_sub(py::module &m) {
    m.def("gbm", &gbm, "Pricing using Stochastic process: GBM");
    m.def("heston", &heston, "Pricing using Stochastic process: Heston");
}


void init_exotic_mont_sub(py::module &m) {
    py::module asian_sub = m.def_submodule("asian", "A nested submodule for asian options");

    init_asian_mont_sub(asian_sub);
}

void init_european_mont_sub(py::module &m) {
    py::module vanilla_sub = m.def_submodule("vanilla", "A nested submodule for vanilla options");
    py::module exotic_sub = m.def_submodule("exotic", "A nested submodule for exotic options");

    init_vanilla_mont_sub(vanilla_sub);
    init_exotic_mont_sub(exotic_sub);

}

void init_montecarlo_sub(py::module &m) {
    py::module european = m.def_submodule("european", "A nested submodule for european style options");
    init_european_mont_sub(european);
}



// submodules for COS method


void init_asian_cos_sub(py::module &m) {
    m.def("gbm", &asian_gbm_1, "Pricing using Stochastic process: GBM");
}


void init_vanilla_cos_sub(py::module &m) {
    m.def("gbm", &gbm_1, "Pricing using Stochastic process: GBM");
    m.def("heston", &heston_1, "Pricing using Stochastic process: Heston");
}


void init_exotic_cos_sub(py::module &m) {
    py::module asian_sub = m.def_submodule("asian", "A nested submodule for asian options");

    init_asian_cos_sub(asian_sub);
}

void init_european_cos_sub(py::module &m) {
    py::module vanilla_sub = m.def_submodule("vanilla", "A nested submodule for vanilla options");
    py::module exotic_sub = m.def_submodule("exotic", "A nested submodule for exotic options");

    init_vanilla_cos_sub(vanilla_sub);
    init_exotic_cos_sub(exotic_sub);

}

void init_cos_method_sub(py::module &m) {
    py::module european = m.def_submodule("european", "A nested submodule for european style options");
    init_european_cos_sub(european);
}

