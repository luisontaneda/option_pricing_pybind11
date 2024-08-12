#include "../COS_method_wrapper.h"
#include "cos_bermudan.cpp"

using namespace std;

namespace COS_method {
namespace american {

complex<double> c_f_gbm(double u, double sigma, double r, double dt) {
  double variance = pow(sigma, 2);
  complex<double> drift = (r - 0.5 * variance) * i_c * u * dt;
  return exp(drift - 0.5 * variance * pow(u, 2) * dt);
}

pair<double, double> a_b_gbm(double x_0, double r, double sigma, double dt,
                             int L) {
  double c1 = (r - .5 * pow(sigma, 2)) * dt;
  double c2 = pow(sigma, 2) * dt;
  double c4 = 0;

  double a = (c1 + x_0) - L * sqrt(c2 + sqrt(c4));
  double b = (c1 + x_0) + L * sqrt(c2 + sqrt(c4));

  return make_pair(a, b);
}

double gbm(double S_0, double K, double T, double r, double sigma, int n_steps,
           int L, int N, int d, string CP) {
  auto c_f = [sigma, r](double u, double dt) -> complex<double> {
    return c_f_gbm(u, sigma, r, dt);
  };

  double x_0 = log(S_0 / K);

  auto a_b_f = [sigma, r, x_0, L](double dt) -> pair<double, double> {
    return a_b_gbm(x_0, r, sigma, dt, L);
  };

  // Extrapolation from bermudan to american option
  double temp_1 =
      cos_bermudan(c_f, a_b_f, S_0, K, T, r, sigma, pow(2, d), L, N, CP);
  double temp_2 =
      cos_bermudan(c_f, a_b_f, S_0, K, T, r, sigma, pow(2, d + 1), L, N, CP);
  double temp_3 =
      cos_bermudan(c_f, a_b_f, S_0, K, T, r, sigma, pow(2, d + 2), L, N, CP);
  double temp_4 =
      cos_bermudan(c_f, a_b_f, S_0, K, T, r, sigma, pow(2, d + 3), L, N, CP);

  double res = (64 * temp_4 - 56 * temp_3 + 14 * temp_2 - temp_1) / 21;

  return res;
}
}  // namespace american
}  // namespace COS_method
