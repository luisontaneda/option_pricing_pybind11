#include <Eigen/Dense>
#include <vector>

//#include "../../helper_functions/h_f_1.cpp"
#include "cos_method.cpp"

using namespace std;

Eigen::MatrixXcd Clenshaw_Curtis_quadrature(double a, double b, int n_q,
                                            int N) {
  int len_mat = n_q / 2 + 1;
  Eigen::MatrixXd D(len_mat, len_mat);
  double temp;

  for (int k = 0; k < len_mat; k++) {
    for (int n = 0; n < len_mat; n++) {
      temp = 2.0 / n_q * cos((n * k * PI) / (n_q / 2.0));

      if (n == 0 || n == len_mat - 1) temp *= .5;

      D(k, n) = temp;
    }
  }

  Eigen::VectorXd d(len_mat);

  d(0) = 1;

  for (int k = 1; k < len_mat - 1; k++) {
    d(k) = 2 / (1 - pow(2 * k, 2));
  }

  d(len_mat - 1) = 1 / (1 - pow(n_q, 2));

  Eigen::MatrixXd w = (D.transpose() * d).transpose();

  auto f_x = [b, a](double x, double u_k, double u_l) -> complex<double> {
    double temp_exp = (b - a) / 2.0 * x + (a + b) / 2.0;
    return (b - a) / 2.0 *
           (pow(exp(temp_exp) + 1, i_c * u_k) * cos((temp_exp - a) * u_l));
  };

  Eigen::VectorXcd y_mat(len_mat);
  double u_k;
  double u_l;
  Eigen::MatrixXcd M_c_c(N, N);

  for (int k = 0; k < N; k++) {
    u_k = k * PI / (b - a);
    for (int l = 0; l < N; l++) {
      u_l = l * PI / (b - a);
      for (int n = 0; n < len_mat; n++) {
        temp = cos(n * PI / n_q);

        y_mat(n) = f_x(temp, u_k, u_l) + f_x(-temp, u_k, u_l);
      }

      M_c_c(k, l) = (w * y_mat)[0];
    }
  }

  return M_c_c;
}

double asian_cos_method(function<complex<double>(double)> c_f, double S_0,
                        double K, double T, double r, int n_steps, int N,
                        int n_q, int L, double a, double b,
                        string CP = "Call") {
  int len_mat = n_q / 2 + 1;
  double dt = T / n_steps;

  Eigen::MatrixXcd M_c_c = Clenshaw_Curtis_quadrature(a, b, n_q, N);

  double u_k;
  complex<double> cf_y[N];
  complex<double> cf_R[N];

  complex<double> temp_cf;

  for (int k = 0; k < N; k++) {
    // first iteration char function
    u_k = k * PI / (b - a);

    temp_cf = c_f(u_k);

    cf_R[k] = temp_cf;
    cf_y[k] = temp_cf;
  }

  Eigen::VectorXcd A(N);
  Eigen::VectorXcd cf_z(N);
  double u_l;

  for (int j = 1; j < n_steps; j++) {
    for (int l = 0; l < N; l++) {
      u_l = l * PI / (b - a);
      A(l) = 2.0 / (b - a) * (cf_y[l] * exp(-i_c * u_l * a)).real();
    };

    A(0) *= .5;
    cf_z = M_c_c * A;

    for (int k = 0; k < N; k++) {
      cf_y[k] = cf_R[k] * cf_z(k);
    }
  }

  double Opt_val;
  double h_k;
  double c;
  double d;

  for (int k = 0; k < N; k++) {
    u_k = k * PI / (b - a);

    c = a;
    d = log(K * (n_steps + 1) / S_0 - 1);

    pair<double, double> chi_psi_pair = chi_psi(a, b, c, d, k);

    h_k = 2.0 / (b - a) *
          ((K - S_0 / (n_steps + 1)) * chi_psi_pair.second -
           (S_0 / (n_steps + 1)) * chi_psi_pair.first);

    if (k == 0) h_k *= .5;

    Opt_val += (cf_y[k] * exp(-i_c * u_k * a)).real() * h_k;
  };

  Opt_val *= exp(-r * T);

  if (CP == "Put") {
    return Opt_val;
  } else {
    double put_call_par = Opt_val +
                          S_0 * exp(-r * T) / (n_steps + 1) *
                              (exp(r * dt * (n_steps + 1)) - 1) /
                              (exp(r * dt) - 1) -
                          K * exp(-r * T);

    return put_call_par;
  }
}