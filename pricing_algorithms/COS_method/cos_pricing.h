#include "asian_cos_pricing.h"

using namespace std;


complex<double> c_f_gbm(double u, double sigma, double r, double tau) {

    double variance = pow(sigma, 2);
    complex<double> drift = (r - 0.5 * variance) * i_c * u * tau;
    return exp(drift - 0.5 * variance * pow(u, 2) * tau);

}


complex<double> c_f_heston(double u, double gmma, double r, double tau, double rho, double v \
                            , double kappa, double v_mac){


    // v_mac: macron v is the expected volatility when t approaches infinity
    // gmma: volatility wich determines variance of v_t
    // 

    complex<double> temp = kappa - i_c*rho*gmma*u;
    complex<double> d_1 = sqrt(pow(temp, 2) + (pow(u, 2) + i_c * u)*pow(gmma, 2));
    complex<double> g = (temp - d_1) / (temp + d_1);

    complex<double> temp_1 = i_c*u*r*tau + (v/pow(gmma, 2)) * (1.0 - exp(-d_1*tau))/(1.0 - g*exp(-d_1*tau)) \
                            * (temp - d_1);

    complex<double> temp_2 = ((kappa*v_mac) / pow(gmma, 2)) * (tau*(temp - d_1) - \
                            2.0*log((1.0 - g*exp(-d_1*tau))/(1.0 - g)));

    return exp(temp_1 + temp_2);

}





complex<double> c_f_ret_gbm(double u, double sigma, double r, double dt) {
    
    // Characteristic Function of the return of a levy procces using
    // gbm as the stochastic process

   return  exp((i_c*u*(r-.5*pow(sigma, 2)) - .5*pow(sigma, 2)*pow(u, 2))*dt);
}


pair<double, double> a_b_interval(double r, double sigma, int time_steps, int L){

    double c_1 = (r - 0.5 * pow(sigma, 2));
    double c_2 = pow(sigma, 2);
    double c_4 = 0;

    double a = 999999;
    double b = -999999;
    double temp_a;
    double temp_b;


    for (int j=1; j<time_steps+1; j++){

        temp_a = c_1 + log(j) - L*sqrt(j*c_2 + sqrt(j*c_4));
        temp_b = j*c_1 + log(j) + L*sqrt(j*c_2 + sqrt(j*c_4));

        if (temp_a < a)
            a = temp_a;

        if (temp_b > b)
            b = temp_b;

    }

    return make_pair(a, b);

}


double gbm_1(double S_0,double K, double tau,double r, double sigma, int L, int N, string CP = "Call"){

    auto   c_f	= [sigma, r, tau] (double u) -> complex<double> {return c_f_gbm(u, sigma, r, tau);};

    return cos_method(c_f, S_0, K, tau, r, L, N, CP);

}



double heston_1(double S_0,double K, double tau, double r, double gmma, double kappa, \
                double rho, double v, double v_mac, int L, int N, string CP = "Call"){


    auto   c_f	= [gmma, r, tau, rho, v, kappa, v_mac] (double u) -> complex<double> {return c_f_heston(u, gmma, r, tau, rho, v \
                                                                                            , kappa, v_mac);};

    return cos_method(c_f, S_0, K, tau, r, L, N, CP);

}


double asian_gbm_1(double S_0,double K, double tau, double r, double sigma, int time_steps, int N, int n_q, int L,string CP = "Call"){

    double dt = tau/time_steps;

    pair<double, double> a_b = a_b_interval(r, sigma, time_steps, L);

    auto   c_f	= [sigma, r, dt] (double u) -> complex<double> {return c_f_ret_gbm(u, sigma, r, dt);};

    return asian_cos_method(c_f, S_0, K, tau, r, time_steps, N, n_q, L, a_b.first, a_b.second, CP);

}