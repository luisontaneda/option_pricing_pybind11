#include <iostream>
#include <random>
#include <string>

#include "montecarlo_wrapper.h"

using namespace std;

int main() {
  int L = 10;

  int N = pow(2, 10);

  // cout << "hola" << endl;

  double S_t = 36.0;
  double K = 40.0;
  double sigma = .4;
  double r = .06;
  double T = 1.0;
  double x_0 = log(S_t / K);
  // int M = 50;
  int n_b_fun = 3;
  int d = 3;
  int n_steps = 100;
  int n_sims = 250000;
  string CP = "Put";

  double res = montecarlo::american::gbm(S_t, K, T, r, sigma, n_b_fun, n_steps, n_sims, CP);

  //double res =
   //   montecarlo::european::gbm(S_t, K, T, r, sigma, n_steps, n_sims, CP);

  
  //double res =
  //    montecarlo::european::gbm(S_t, K, T, r, sigma, n_steps, n_sims, CP);

  return 0;
}

