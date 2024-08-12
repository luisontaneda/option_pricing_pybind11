#ifndef COS_METHOD_FUNCTIONS
#define COS_METHOD_FUNCTIONS

#include <iostream>
#include <string>

namespace COS_method {
namespace american {

double gbm(double S_0, double K, double T, double r, double sigma,
           int n_steps, int L, int N, int d = 3, std::string CP = "Call");
           
}  // namespace american

namespace european {

double gbm(double S_0, double K, double T, double r, double sigma, int L,
           int N, std::string CP = "Call");

double heston(double S_0, double K, double T, double r, double gmma,
              double kappa, double rho, double v, double v_mac, int L, int N,
              std::string CP = "Call");

double asian_gbm(double S_0, double K, double T, double r, double sigma,
                 int n_steps, int N, int n_q, int L,
                 std::string CP = "Call");

}  // namespace european
}  // namespace COS_method

#endif  // MONTECARLO_FUNCTIONS
