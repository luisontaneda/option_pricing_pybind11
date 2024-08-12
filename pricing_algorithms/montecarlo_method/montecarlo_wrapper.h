#ifndef MONTECARLO_FUNCTIONS
#define MONTECARLO_FUNCTIONS

#include <iostream>
#include <string>

namespace montecarlo {
namespace american {

double gbm(double S_t, double K, double T, double r, double sigma,
           int n_b_fun = 3, int n_steps = 100, int n_sims = 250000,
           std::string CP = "Call");

}  // namespace american

namespace european {

double gbm(double S_t, double K, double T, double r, double sigma,
           int n_steps = 100, int n_sims = 250000, std::string CP = "Call");

double heston(double S_t, double K, double T, double r, double gmma,
              double kappa, double rho, double v_i, double v_mac,
              int n_steps = 100, int n_sims = 250000, std::string CP = "Call");

double asian_gbm(double S_t, double K, double T, double r, double sigma,
                 int n_steps = 100, int n_sims = 250000,
                 std::string CP = "Call");

}  // namespace european
}  // namespace montecarlo

#endif  // MONTECARLO_FUNCTIONS
