#include <Eigen/Dense>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>

// #include "asian_cos_pricing.h"
#include "../../helper_functions/h_f_1.h"

using namespace std;

Eigen::VectorXcd circularConvolution(const Eigen::VectorXcd& m,
                                     const Eigen::VectorXcd& u) {
  int N = m.size();
  Eigen::VectorXcd result(N);

  for (int n = 0; n < N; ++n) {
    result(n) = 0;
    for (int k = 0; k < N; ++k) {
      int index = (n - k + N) % N;
      result(n) += m(k) * u(index);
    }
  }

  return result;
}

void func_FFT(function<complex<double>(double)> c_f, double* V_k,
              Eigen::VectorXcd& M_s_u, Eigen::VectorXcd& M_c_u, double x_2,
              float x_1, double a, double b, int N) {
  Eigen::VectorXcd m_s(2 * N);
  Eigen::VectorXcd m_c(2 * N);
  Eigen::VectorXcd u_s(2 * N);
  Eigen::VectorXcd u_c(2 * N);
  double j_double;
  double u_temp;
  int idx_j;

  m_s(0) = i_c * ((x_2 - x_1) * PI) / (b - a);
  j_double = 2 * N - 1;
  m_c(0) = (exp(i_c * j_double * (x_2 - a) * PI / (b - a)) -
            exp(i_c * j_double * (x_1 - a) * PI / (b - a))) /
           j_double;

  u_s(0) = .5 * c_f(0) * V_k[0];
  u_c(N) = .5 * c_f(0) * V_k[0];

  for (int j = 1; j < 2 * N; j++) {
    j_double = 2 * N - (j + 1);
    m_c(j) = (exp(i_c * j_double * (x_2 - a) * PI / (b - a)) -
              exp(i_c * j_double * (x_1 - a) * PI / (b - a))) /
             j_double;

    if (j < N) {
      j_double = (2 * N - j) - 2 * N;
      u_temp = j * PI / (b - a);
      u_s(j) = c_f(u_temp) * V_k[j];
    } else if (j > N) {
      idx_j = j - N;
      u_temp = idx_j * PI / (b - a);
      u_c(j) = c_f(u_temp) * V_k[idx_j];

      j_double = 2 * N - j;
    }

    m_s(j) = (exp(i_c * j_double * (x_2 - a) * PI / (b - a)) -
              exp(i_c * j_double * (x_1 - a) * PI / (b - a))) /
             j_double;
  }

  m_c(2 * N - 1) = i_c * ((x_2 - x_1) * PI) / (b - a);
  m_s(N) = 0;

  M_s_u = circularConvolution(m_s, u_s).head(N);
  M_c_u = circularConvolution(m_c, u_c).head(N).reverse();
}

double get_root(string CP, double guess, function<double(double)> func,
                function<double(double)> f_prime,
                uniform_real_distribution<> dis, double low_b, double up_b) {
  random_device rd;
  mt19937 gen(rd());

  for (int i = 0; i < 50; i++) {
    try {
      if (i > 0) {
        if (CP == "Put") {
          guess = dis(gen);
        } else {
          guess = dis(gen);
        }
      }

      double x_ast = boost::math::tools::newton_raphson_iterate(
          [&](double x) {
            return std::make_pair(func(x), f_prime(x));
          },                                   // Function and derivative
          guess,                               // Initial guess
          low_b,                               // Lower bound for the root
          up_b,                                // Upper bound for the root
          std::numeric_limits<double>::digits  // Maximum precision
      );
      return x_ast;
    } catch (const std::exception& e) {
      std::cerr << "Error with initial guess " << guess << ": " << e.what()
                << std::endl;
    }
  }
  return guess;
}

double cos_bermudan(function<complex<double>(double, double)> c_f_1,
                    function<pair<double, double>(double)> a_b_f, double S_0,
                    double K, double T, double r, double sigma, int n_steps,
                    int L, int N, string CP = "Call") {
  double x_0 = log(S_0 / K);
  double dt = T / n_steps;

  pair<double, double> a_b_pair = a_b_f(dt);
  double a = a_b_pair.first;
  double b = a_b_pair.second;

  auto c_f = [c_f_1, dt](double u) -> complex<double> { return c_f_1(u, dt); };
  auto G = [a, b, K, CP](double k, double x_1, double x_2) -> double {
    double H_k;

    if (CP == "Put") {
      double c = x_1;
      double d = x_2;
      pair<double, double> chi_psi_pair = chi_psi(a, b, c, d, k);
      H_k = 2.0 / (b - a) * (-chi_psi_pair.first + chi_psi_pair.second);
    } else {
      double c = x_1;
      double d = x_2;
      pair<double, double> chi_psi_pair = chi_psi(a, b, c, d, k);
      H_k = 2.0 / (b - a) * (chi_psi_pair.first - chi_psi_pair.second);
    }
    return K * H_k;
  };

  double V_t[N] = {0};

  double alpha_CP;

  if (CP == "Put") {
    for (int k = 0; k < N; k++) {
      V_t[k] = G(k, a, 0);
    }

    alpha_CP = -1;
  } else {
    for (int k = 0; k < N; k++) {
      V_t[k] = G(k, 0, b);
    }

    alpha_CP = 1;
  }

  Eigen::VectorXcd M_s_u;
  Eigen::VectorXcd M_c_u;

  auto f_prime = [r, dt, K, N, a, b, alpha_CP, c_f, &V_t](double x) -> double {
    double temp_1 = 0;
    double u_temp;

    for (int k = 0; k < N; k++) {
      u_temp = k * PI / (b - a);
      temp_1 +=
          (i_c * u_temp * c_f(u_temp) * exp(i_c * u_temp * (x - a))).real() *
          V_t[k];

      if (k == 0) temp_1 *= .5;
    }

    temp_1 *= exp(-r * dt);
    temp_1 -= alpha_CP * K * exp(x);

    return temp_1;
  };

  auto func = [r, dt, K, N, a, b, alpha_CP, c_f, &V_t](double x) -> double {
    double temp_1 = 0;
    double u_temp;

    for (int k = 0; k < N; k++) {
      u_temp = k * PI / (b - a);
      temp_1 += (c_f(u_temp) * exp(i_c * u_temp * (x - a))).real() * V_t[k];

      if (k == 0) temp_1 *= .5;
    }

    temp_1 *= exp(-r * dt);
    temp_1 -= alpha_CP * K * (exp(x) - 1);

    return temp_1;
  };

  Eigen::MatrixXd C_approx;

  double x_ast = (a + b) / 2;

  for (int m = n_steps - 1; m > 0; m--) {
    uniform_real_distribution<> dis;

    if (CP == "Put") {
      dis.param(uniform_real_distribution<>::param_type(x_ast, 0.0));
      x_ast = get_root(CP, x_ast, func, f_prime, dis, x_ast - 10, 0);
      if (x_ast < a) {
        x_ast = a;
      }
    } else {
      dis.param(uniform_real_distribution<>::param_type(0.0, x_ast));
      x_ast = get_root(CP, x_ast, func, f_prime, dis, 0, x_ast + 10);
      if (x_ast > b) {
        x_ast = b;
      }
    }

    if (CP == "Put") {
      func_FFT(c_f, V_t, M_s_u, M_c_u, b, x_ast, a, b, N);
      C_approx = exp(-r * dt) * (M_s_u + M_c_u).imag() / PI;

      for (int k = 0; k < N; k++) {
        V_t[k] = C_approx(k) + G(k, a, x_ast);
      }
    } else {
      func_FFT(c_f, V_t, M_s_u, M_c_u, x_ast, a, a, b, N);
      C_approx = exp(-r * dt) * (M_s_u + M_c_u).imag() / PI;

      for (int k = 0; k < N; k++) {
        V_t[k] = C_approx(k) + G(k, x_ast, b);
      }
    }
  }

  double u_temp;
  double V_fin = 0;

  for (int k = 0; k < N; k++) {
    u_temp = k * PI / (b - a);
    V_fin += (c_f(u_temp) * exp(i_c * u_temp * (x_0 - a))).real() * V_t[k];

    if (k == 0) V_fin *= .5;
  }

  V_fin *= exp(-r * dt);

  return V_fin;
}