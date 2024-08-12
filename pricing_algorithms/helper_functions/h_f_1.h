//#include <functional>
#include <cmath>
#include <complex>
#include <iostream>

using namespace std;

const complex<double> i_c(0.0, 1.0);
const double PI = 3.14159265358979323846;

static pair<double, double> chi_psi(double a, double b, double c, double d, double k) {
  double psi =
      sin(k * PI * (d - a) / (b - a)) - sin(k * PI * (c - a) / (b - a));

  if (k == 0)
    psi = d - c;
  else
    psi = psi * (b - a) / (k * PI);

  double chi = 1.0 / (1.0 + pow((k * PI / (b - a)), 2));

  double expr1 = cos(k * PI * (d - a) / (b - a)) * exp(d) -
                 cos(k * PI * (c - a) / (b - a)) * exp(c);

  double expr2 = k * PI / (b - a) * sin(k * PI * (d - a) / (b - a)) * exp(d) -
                 k * PI / (b - a) * sin(k * PI * (c - a) / (b - a)) * exp(c);
  chi = chi * (expr1 + expr2);

  return make_pair(chi, psi);
}

static double H_k_func(double a, double b, double k, string CP = "Call") {
  double c;
  double d;
  double H_k;

  if (CP == "Call") {
    c = 0.0;
    d = b;

    pair<double, double> chi_psi_pair = chi_psi(a, b, c, d, k);

    H_k = 2.0 / (b - a) * (chi_psi_pair.first - chi_psi_pair.second);
  } else if (CP == "Put") {
    c = a;
    d = 0.0;

    pair<double, double> chi_psi_pair = chi_psi(a, b, c, d, k);

    H_k = 2.0 / (b - a) * (-chi_psi_pair.first + chi_psi_pair.second);
  }

  return H_k;
}

