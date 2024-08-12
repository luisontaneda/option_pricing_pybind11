#include "cos_american.cpp"

int main() {
  int L = 10;

  int N = pow(2, 11);

  double S_t = 36.0;
  double K = 40.0;
  double sigma = .4;
  double r = .06;
  double T = 1.0;
  double x_0 = log(S_t / K);
  string CP = "Call";
  // int M = 50;
  int d = 2;

  double temp_1 = cos_bermudan(S_t, K, T, r, sigma, pow(2, d), L, N, CP);
  double temp_2 = cos_bermudan(S_t, K, T, r, sigma, pow(2, d + 1), L, N, CP);
  double temp_3 = cos_bermudan(S_t, K, T, r, sigma, pow(2, d + 2), L, N, CP);
  double temp_4 = cos_bermudan(S_t, K, T, r, sigma, pow(2, d + 3), L, N, CP);

  double res = (64 * temp_4 - 56 * temp_3 + 14 * temp_2 - temp_1) / 21;

  return 0;
}