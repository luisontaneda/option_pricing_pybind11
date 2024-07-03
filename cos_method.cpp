#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#include <functional>

using namespace std;

const complex<double> i_c(0.0,1.0);
const double PI = 3.14159265358979323846;

pair<double, double> chi_psi(double a, double b, double c, double d, double k){

    double psi = sin(k * PI * (d - a) / (b - a)) - sin(k * PI * (c - a)/(b - a));
    
    if (k == 0)
        psi = d - c;
    else
        psi = psi * (b - a) / (k *  PI);


    double chi = 1.0 / (1.0 + pow((k * PI / (b - a)) , 2));

    double expr1 = cos(k * PI * (d - a)/(b - a)) * exp(d)  - cos(k * PI \
                  * (c - a) / (b - a)) * exp(c);

    double expr2 = k * PI / (b - a) * sin(k * PI * 
                        (d - a) / (b - a)) * exp(d)  - k * PI / (b - a) * sin(k \
                        * PI * (c - a) / (b - a)) * exp(c);
    chi = chi * (expr1 + expr2);

    return make_pair(chi, psi);

}

double H_k_func(double a, double b, double k, string CP = "Call"){

    double c;
    double d;
    double H_k;

    if (CP == "Call"){

        c = 0.0;
        d = b;

        pair<double, double> chi_psi_pair = chi_psi(a, b, c, d, k);

        H_k = 2.0 / (b - a) * (chi_psi_pair.first - chi_psi_pair.second);
    }
    else if (CP == "Put"){

        c = a;
        d = 0.0;

        pair<double, double> chi_psi_pair = chi_psi(a, b, c, d, k);

        H_k = 2.0 / (b - a) * (-chi_psi_pair.first + chi_psi_pair.second);
    }

    return H_k;
}


double cos_method(function<complex<double>(double)> c_f, double S_0,double K, double tau, double r, int L, int N = pow(2,9), string CP = "Call")
{   

    // L - Size of truncation domain
    // tau - time to maturity
    //


    double x_0 = log(S_0/K);

    // Integration (truncation) domain
    double a = 0 - L * sqrt(tau);
    double b = 0 + L * sqrt(tau);

    double u;
    double h_k;
    double Opt_val;
    complex<double> cf_ev;

    for (int k = 0; k < N; k++){

        u = k * PI/(b-a);
        cf_ev = c_f(u);
        h_k = H_k_func(a, b, k, CP);

        if (k == 0)
            h_k *= .5;

        Opt_val += (cf_ev * exp(u * i_c * (x_0 - a))).real() * h_k;

    };

    Opt_val *= exp(-r * tau) * K;
    
    return Opt_val;
}