#include <pybind11/pybind11.h>
#include "pricing_algorithms/COS_method/COS_method_wrapper.h"
#include "pricing_algorithms/montecarlo_method/montecarlo_wrapper.h"

namespace py = pybind11;

void init_montecarlo(py::module&);
void init_cos_method(py::module&);

PYBIND11_MODULE(option_pricing, m) {
    m.doc() = "option_pricing plugin";

    py::module montecarlo = m.def_submodule("montecarlo", "Submodule for Monte Carlo methods");
    py::module cos_method = m.def_submodule("cos_method", "Submodule for COS method");

    init_montecarlo(montecarlo);
    init_cos_method(cos_method);
}

// Monte Carlo method submodules

void init_montecarlo(py::module& m) {
    py::module european = m.def_submodule("european", "Submodule for European style options");
    {
    py::module vanilla = european.def_submodule("vanilla", "Submodule for Vanilla options");
    py::module exotic = european.def_submodule("exotic", "Submodule for Exotic options");
    py::module asian = exotic.def_submodule("asian", "Submodule for Asian options");

    vanilla.def("gbm", &montecarlo::european::gbm, "Pricing using Stochastic process: GBM");
    vanilla.def("heston", &montecarlo::european::heston, "Pricing using Stochastic process: Heston");
    asian.def("gbm", &montecarlo::european::asian_gbm, "Pricing using Stochastic process: GBM");
    }

    py::module american = m.def_submodule("american", "Submodule for American style options");
    {
    py::module vanilla = american.def_submodule("vanilla", "Submodule for Vanilla options");

    vanilla.def("gbm", &montecarlo::american::gbm, "Pricing using Stochastic process: GBM");
    }
}

// COS method submodules

void init_cos_method(py::module& m) {
    py::module european = m.def_submodule("european", "Submodule for European style options");
    {
    py::module vanilla = european.def_submodule("vanilla", "Submodule for Vanilla options");
    py::module exotic = european.def_submodule("exotic", "Submodule for Exotic options");
    py::module asian = exotic.def_submodule("asian", "Submodule for Asian options");

    vanilla.def("gbm", &COS_method::european::gbm, "Pricing using Stochastic process: GBM");
    vanilla.def("heston", &COS_method::european::heston, "Pricing using Stochastic process: Heston");
    asian.def("gbm", &COS_method::european::asian_gbm, "Pricing using Stochastic process: GBM");
    }

    py::module american = m.def_submodule("american", "Submodule for American style options");
    {
    py::module vanilla = american.def_submodule("vanilla", "Submodule for Vanilla options");

    vanilla.def("gbm", &COS_method::american::gbm, "Pricing using Stochastic process: GBM");
    }
}