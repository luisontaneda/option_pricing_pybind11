#include <string>
#include "montecarlo_pricing.h"


double gbm_process(mt19937 gen, double S_t, double r, double sigma, double dt, double n_steps){


    normal_distribution<double> Z(0.0, 1.0);
    double x_i = log(S_t);

    for (int j=0;j<=n_steps;j++)
    {
        x_i += (r - (pow(sigma, 2)/2))*dt + sigma * sqrt(dt) * Z(gen);
    }

    return exp(x_i);
}


double heston_process(mt19937 gen \
                     ,double S_t, double r, double dt, double kappa, double gmma, \
                        double v_mac, double rho, double v_i, int n_steps){

    
    double x_i = log(S_t);
    double v_i_1;
    double c_mac = (pow(gmma, 2)/(4*kappa)) * (1 - exp(-kappa*(dt)));
    double delta = (4*kappa*v_mac)/pow(gmma, 2);
    double k_0 = (r - (rho/gmma)*kappa*v_mac)*dt;
    double k_1 = ((rho*kappa)/gmma - .5)*dt - (rho/gmma);
    double k_2 = rho/gmma;
    double k_3 = (1 - pow(rho, 2))*dt;
    double X_2;

    auto k_mac	= [kappa, dt, gmma] (double v_i) -> double {return (4*kappa*exp(-kappa*dt)*v_i)/(pow(gmma, 2)*(1 - exp(-kappa*dt)));};

    normal_distribution<double> Z(0.0, 1.0);
    
    double df;


    for (int j=0;j<=n_steps;j++){

        df = k_mac(v_i);

        poisson_distribution<int> pois(df/2);
        X_2 = pois(gen);

        if (X_2 > 0){
            chi_squared_distribution<double> chi_2(2*X_2);
            X_2 = chi_2(gen);
        }

        if (df > 0){
            chi_squared_distribution<double> chi(delta);
            X_2 += chi(gen);
        }

        v_i_1 = c_mac*X_2;
        
        x_i += k_0 + k_1*v_i + k_2*v_i_1 + sqrt(k_3*v_i)*Z(gen);

        v_i = v_i_1;

    };

    return exp(x_i);
}



double gbm_process_1(mt19937 gen, double S_t, double r, double sigma, double dt, int n_steps){


    normal_distribution<double> Z(0.0, 1.0);
    double arr[n_steps+1];
    arr[0] = S_t;

    for (int j=1;j<=n_steps+1;j++)
    {
        S_t = S_t*exp((r - pow(sigma, 2)/2)*dt + sigma * sqrt(dt) * Z(gen));
        arr[j] = S_t;
    }


    return calculate_mean(arr, n_steps+1);
}





double gbm(double S_t, double K, double tau, double r, double sigma, int n_steps = 100, int n_sims = 250000, string CP = "Call"){


    double dt = tau/n_steps;

    auto   stoch_p	= [S_t, r, dt, sigma, n_steps] (mt19937 gen) -> double {return gbm_process(gen, S_t, r, sigma, dt, n_steps);};

    return m_c_pricing(stoch_p, S_t, K, tau, r, n_sims, CP);    
}



double heston(double S_t, double K, double tau, double r, double gmma, double kappa \
            , double rho, double v_i, double v_mac, int n_steps = 100, int n_sims = 250000, string CP = "Call"){


    double dt = tau/n_steps;

    auto   stoch_p	= [S_t, r, dt, kappa, gmma, v_mac, rho, v_i, n_steps] (mt19937 gen) -> double {return heston_process(gen, S_t, r, dt, kappa, gmma, v_mac, rho, v_i, n_steps);};

    return m_c_pricing(stoch_p, S_t, K, tau, r, n_sims, CP);    
}



double asian_gbm(double S_t, double K, double tau, double r, double sigma, int n_steps = 100, int n_sims = 250000,string CP = "Call"){

    double dt = tau/n_steps;

    auto   stoch_p	= [S_t, r, dt, sigma, n_steps] (mt19937 gen) -> double {return gbm_process_1(gen, S_t, r, sigma, dt, n_steps);};

    return m_c_pricing(stoch_p, S_t, K, tau, r, n_sims, CP);

}
