#include <iostream>
#include <random>
#include <string>

#include "COS_method_wrapper.h"

using namespace std;

int main() {
  int L = 10;

  int N = pow(2, 9);

  // cout << "hola" << endl;

  double S_0 = 36.0;
  double K = 40.0;
  double sigma = .4;
  double r = .06;
  double T = 1.0;
  double x_0 = log(S_0 / K);
  // int M = 50;
  int n_b_fun = 3;
  int d = 3;
  int n_steps = 100;
  int n_sims = 250000;
  int n_q = 1000;
  string CP = "Put";

  double res =
      COS_method::american::gbm(S_0, K, T, r, sigma, n_steps, L, N, d, CP);

  // double res = COS_method::european::gbm(S_0, K, T, r, sigma, L, N, CP);

  // double res = COS_method::european::asian_gbm(S_0, K, T, r, sigma,
  //                n_steps, N, n_q, L, CP);

  return 0;
}