#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <string>

#include "../../helper_functions/h_f_1.h"

using namespace std;

double cos_method(function<complex<double>(double)> c_f, double S_0, double K,
                  double T, double r, int L, int N = pow(2, 9),
                  string CP = "Call") {

  // S_0 - Initial price
  // K - Strike price
  // T - Time to maturity
  // r - Expected rate of return
  // L - Size of truncation domain
  // N - Number of cosine terms in the expansion
  // CP - Option type

  double x_0 = log(S_0 / K);

  // Integration (truncation) domain
  double a = 0 - L * sqrt(T);
  double b = 0 + L * sqrt(T);

  double u;
  double h_k;
  double Opt_val;
  complex<double> cf_ev;

  for (int k = 0; k < N; k++) {
    u = k * PI / (b - a);
    cf_ev = c_f(u);
    h_k = H_k_func(a, b, k, CP);

    if (k == 0) h_k *= .5;

    Opt_val += (cf_ev * exp(u * i_c * (x_0 - a))).real() * h_k;
  };

  Opt_val *= exp(-r * T) * K;

  return Opt_val;
}