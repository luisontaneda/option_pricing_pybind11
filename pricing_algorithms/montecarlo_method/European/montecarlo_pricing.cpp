#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

using namespace std;

double calculate_mean(const double* numbers, double size) {
  if (size == 0) return 0.0f;  // Handle empty array case

  double sum = 0.0f;
  for (int i = 0; i < size; ++i) {
    sum += numbers[i];
  }
  return sum / size;
}

double m_c_pricing(function<double(mt19937)> stoch_p, double S_t, double K,
                   double T, double r, int n_sims, string CP) {
  // S_t: Stock price at day 0
  // K: Strike price
  // T: Expiry time horizon in years
  //

  double arr[n_sims];

  for (int i = 0; i < n_sims; i++) {
    random_device rd;
    mt19937 gen{rd()};

    S_t = stoch_p(gen);

    if (CP == "Call") {
      arr[i] = max<double>(S_t - K, 0.0);
    } else if (CP == "Put") {
      arr[i] = max<double>(-S_t + K, 0.0);
    }
  };

  double mean_sims = exp(-r * T) * calculate_mean(arr, n_sims);

  return mean_sims;
}