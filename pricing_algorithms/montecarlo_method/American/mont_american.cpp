#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <string>

using namespace std;

Eigen::VectorXd la_guerre_pol(double X, int n_b_fun) {
  Eigen::VectorXd L(n_b_fun);

  L(0) = 1;
  L(1) = 1 - X;

  for (int k = 1; k < n_b_fun - 1; ++k) {
    L(k + 1) = ((2 * k + 1 - X) * L(k) - k * L(k - 1)) / (k + 1);
  }

  return L;
}

double m_c_a_pricing(function<int(mt19937, double*)> stoch_p, double S_t,
                     double K, double T, double r, int n_b_fun, int n_steps,
                     int n_sims, string CP) {
  double dt = T / n_steps;

  n_steps += 1;

  double arr_S[n_steps];
  vector<vector<double>> arr_all_S(n_sims, vector<double>(n_steps));
  double final_p_f_prices[n_sims];
  Eigen::VectorXd Y(n_sims);

  double p_f;

  auto func_pay_off = [CP, K](double S) -> double {
    if (CP == "Call") {
      return max<double>(S - K, 0.0);
    } else {
      return max<double>(K - S, 0.0);
    }
  };

  for (int i = 0; i < n_sims; i++) {
    random_device rd;
    mt19937 gen{rd()};

    // the T - dt (second to last) price is used to measure payoff

    stoch_p(gen, arr_S);

    p_f = func_pay_off(arr_S[n_steps - 1]);
    final_p_f_prices[i] = p_f * exp(-r * T);
    Y(i) = p_f;

    for (int k = 0; k < n_steps; k++) {
      arr_all_S[i][k] = arr_S[k];
    };
  };

  Eigen::VectorXd w;
  Eigen::MatrixXd X_mat(n_sims, n_b_fun);
  Eigen::VectorXd cont_exp_X;

  double pay_off_arr[n_sims];
  double temp;

  vector<int> idx_pay_off;

  for (int c = n_steps - 2; c > 0; c--) {
    idx_pay_off.clear();

    for (int i = 0; i < n_sims; i++) {
      temp = arr_all_S[i][c];
      p_f = func_pay_off(temp);
      pay_off_arr[i] = p_f;

      // in the money
      if (p_f > 0) {
        X_mat.row(i) = la_guerre_pol(temp, n_b_fun);
        idx_pay_off.push_back(i);
      }
      // out of the money
      else {
        X_mat.row(i).setZero();
      }
    }

    w = (X_mat.transpose() * X_mat).ldlt().solve(X_mat.transpose() * Y);
    cont_exp_X = X_mat * w;

    for (int i : idx_pay_off) {
      if (cont_exp_X(i) < pay_off_arr[i]) {
        final_p_f_prices[i] = pay_off_arr[i] * exp(-r * dt * c);
        Y(i) = pay_off_arr[i];
      }
    }
  }

  double sum_1;
  for (int i = 0; i < n_sims; i++) {
    sum_1 += final_p_f_prices[i];
  }

  double mean = sum_1 / n_sims;

  return mean;
}