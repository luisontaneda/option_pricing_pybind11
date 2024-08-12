#include "../montecarlo_wrapper.h"
#include "mont_american.cpp"

int gbm_process(mt19937 gen, double S_t, double r, double sigma, double dt,
                int n_steps, double* arr) {
  normal_distribution<double> Z(0.0, 1.0);
  arr[0] = S_t;

  for (int j = 1; j < n_steps+1; j++) {
    S_t = S_t * exp((r - pow(sigma, 2) / 2) * dt + sigma * sqrt(dt) * Z(gen));
    arr[j] = S_t;
  }

  return 0;
}

namespace montecarlo {
namespace american {

double gbm(double S_t, double K, double T, double r, double sigma, int n_b_fun,
           int n_steps, int n_sims, std::string CP) {
  double dt = T / n_steps;

  auto stoch_p = [S_t, r, dt, sigma, n_steps](mt19937 gen,
                                              double* arr_S) -> int {
    return gbm_process(gen, S_t, r, sigma, dt, n_steps, arr_S);
  };

  return m_c_a_pricing(stoch_p, S_t, K, T, r, n_b_fun, n_steps, n_sims, CP);
};
}  // namespace american
}  // namespace montecarlo
