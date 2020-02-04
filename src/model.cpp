#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <string>
#include "brent.hpp"
#include <random> // for distributions
#include <omp.h>
#include <stdio.h>
#include <vector>
#include <algorithm>  
#include <stdexcept> // for handling exceptions
#include "boost/math/special_functions/gamma.hpp" // for the incomplete gamma function
#include "model.hpp"
#include "utils.hpp"

using namespace std;
using namespace brent; // for zero()

using boost::math::tgamma_lower;


// set up random number generators
extern mt19937 generator;

const double check_tolerance(1e-15);

double Teqn(double Ts, double params[]) {
  double A = params[0];
  double b = params[1];
  double C = params[2];
  double n = params[3];
  double s0 = params[4];
  double k = params[5];
  double mu = params[6];

  double As = A * k*Ts;
  double Cs = C * k*Ts;

  double intfac = 1 / mu * pow(Cs, n) * Ts / (n + 1) * pow(1 - s0, n + 1);
  double bnfr = (b + 1) / (n + 1);

  double result = pow(As, b) / pow(Cs, n * bnfr) * pow(mu / Ts * (n + 1), (b - n) / (n + 1)) * tgamma_lower(bnfr, intfac) - exp(-intfac);

  return result;
}

// Ramp load Ts solution - phase2
double Teqn2 ( double Tf, double params[] ) {
  double A = params[0];
  double b = params[1];
  double C = params[2];
  double n = params[3];
  double s0 = params[4];
  double k = params[5];
  double mu = params[6];

  // need to inclue Ts and alpha_T1  
  double Ts = params[7];
  double alpha_T1 = params[8];
  double As = A * k*Ts;
  double Cs = C * k*Ts;
  double intfac  =  1/mu * pow(Cs, n) * Ts/(n+1) * pow(Tf/Ts-s0, n+1);
  double bnfr = (b + 1) / (n + 1);
  
  double result =  pow(As, b) / pow(Cs, n * bnfr) * pow(mu / Ts * (n + 1), (b - n) / (n + 1)) * tgamma_lower(bnfr, intfac) - exp(-intfac) + alpha_T1;
  
  return result;
}

// compute the constants for Tc

vector<double> const_Tc(double Ts, double T0, double params[]) {
  double A = params[0];
  double b = params[1];
  double C = params[2];
  double n = params[3];
  double s0 = params[4];
  double k = params[5];
  double mu = params[6];

  double As = A * k*Ts;
  double Cs = C * k*Ts;

  double intfac = 1 / mu * pow(Cs, n) * Ts / (n + 1) * pow(T0 / Ts - s0, n + 1);
  double bnfr = (b + 1) / (n + 1);

  vector<double> res(4);
  // the constants are C1, C2, alpha_T0 and H^star(T0)

  res[0] = 1 / mu * pow(As * (T0 / Ts - s0), b);
  res[1] = 1 / mu * pow(Cs * (T0 / Ts - s0), n);
  res[2] = exp(intfac) * pow(As, b) / pow(Cs, n * bnfr) * pow(mu / Ts * (n + 1), (b - n) / (n + 1)) * tgamma_lower(bnfr, intfac);
  res[3] = exp(-res[1] * T0);

  return res;
}


double TeqnR ( double Tf, double params[] ) { // ramp load failure time solution
  double A = params[0];
  double b = params[1];
  double C = params[2];
  double n = params[3];
  double s0 = params[4];
  double k_s = params[5];
  double mu = params[6];

  // need to include Ts and applied rate
  double Ts = params[7];
  double k = params[8];
  
  double tau_s = k_s * Ts;
  double As = A * tau_s;
  double Cs = C * tau_s;
  double intfac  =  1/mu * pow(Cs, n) * (tau_s/k)/(n+1) * pow(Tf*k/tau_s-s0, n+1);
  double bnfr = (b + 1) / (n + 1);
  
  double result =  pow(As, b) / pow(Cs, n * bnfr) * pow(mu*k / tau_s * (n + 1), (b - n) / (n + 1)) * tgamma_lower(bnfr, intfac) - exp(-intfac);
  
  return result;
}

// Ramp load with different rate
vector<double> CANR(vector<double> params, int N, double k, double k_s) {  
  double mu_A = params[0];
  double s_A = params[1];
  double mu_b = params[2];
  double s_b = params[3];
  double mu_C = params[4];
  double s_C = params[5];
  double mu_n = params[6];
  double s_n = params[7];
  double mu_s0 = params[8];
  double s_s0 = params[9];

  double t = r8_epsilon();

  double mu_s = 1;  
  
  vector<double> res(N);

  // Generate random variates
  std::normal_distribution<double> normA(mu_A, s_A);
  std::normal_distribution<double> normb(mu_b, s_b);
  normal_distribution<double> normC(mu_C, s_C);
  normal_distribution<double> normn(mu_n, s_n);
  normal_distribution<double> norms0(mu_s0, s_s0);

  omp_set_dynamic(0);
#pragma omp parallel for num_threads(16)
  for (int i = 0; i < N; i++) {
    double tpar[10];
    tpar[0] = exp(normA(generator));
    tpar[1] = exp(normb(generator));
    tpar[2] = exp(normC(generator));
    tpar[3] = exp(normn(generator));
    double tempx = exp(norms0(generator));
    tpar[4] = tempx / (1.0 + tempx);
    tpar[5] = k_s; // the standard loading rate
    tpar[6] = mu_s; //mu_s can be set to an appropriate constant (not estimated).
   
    try {
      double Ts = zero(0.00001, 0.4, t, Teqn, tpar); 
      
      if ((Teqn(0.00001, tpar) > 0) - (Teqn(0.4, tpar) > 0) == 0) {
        res[i] = 1e7; // no solution within usual timeframe (0.1 hours = 38840psi) so set to upbrd.  
        tpar[8] = 0.0;
      } else {
        tpar[7] = Ts;
        tpar[8] = k;

        double Tr = zero(tpar[4]* Ts * k_s/k + 1e-7, 100000/k, t, TeqnR, tpar);
        res[i] = Tr;
        if (TeqnR(Tr, tpar) > 0.1) // numerical error no solution
          res[i] = 1e7;
      
      }
      // manual checking 
//      cout << setprecision(10) << "a=" << tpar[0] << ", b=" << tpar[1] << ", c=" << tpar[2] << ", n=" << tpar[3] << ", s0=" << tpar[4] << " " ;
//      for (int j = 5; j < 9; j++)
//        cout << tpar[j] << " ";
//      cout << res[i] << " " << res[i] * k / (Ts * k_s) << endl;
            
    
    } catch (const exception& e) {
      //			cout << "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
      //			cout << tpar[0] << " " << tpar[1] << " " << tpar[2] << " " << tpar[3] << " " << tpar[4] << " " << tpar[5] << " " << tpar[6] << endl;
    }
  }
  return res;

}

vector<double> CANL(vector<double> params, int N, double tauc, double k, double T1) {
  // tauc is the constant load; if it is negative, it is a ramp load test
  double mu_A = params[0];
  double s_A = params[1];
  double mu_b = params[2];
  double s_b = params[3];
  double mu_C = params[4];
  double s_C = params[5];
  double mu_n = params[6];
  double s_n = params[7];
  double mu_s0 = params[8];
  double s_s0 = params[9];

  double t = r8_epsilon();

  double mu_s = 1;
  double T0 = tauc / k;

  const double uprbd = 1e9; //if the generate number is nan, the result is set to be this number (needed to be increased since 4 years is already larger than 20000)

  double tmp;
  double alpha_T0, alpha_T1;

  vector<double> res(N);

  // Generate random variates
  std::normal_distribution<double> normA(mu_A, s_A);
  std::normal_distribution<double> normb(mu_b, s_b);
  normal_distribution<double> normC(mu_C, s_C);
  normal_distribution<double> normn(mu_n, s_n);
  normal_distribution<double> norms0(mu_s0, s_s0);

  omp_set_dynamic(0);
#pragma omp parallel for num_threads(16)
  for (int i = 0; i < N; i++) {
    double tpar[10];
    tpar[0] = exp(normA(generator));
    tpar[1] = exp(normb(generator));
    tpar[2] = exp(normC(generator));
    tpar[3] = exp(normn(generator));
    double tempx = exp(norms0(generator));
    tpar[4] = tempx / (1.0 + tempx);
    tpar[5] = k; // the loading rate
    tpar[6] = mu_s; //mu_s can be set to an appropriate constant (not estimated).

    try {
      double Ts = zero(0.00001, 0.4, t, Teqn, tpar); 
      tpar[7] = Ts;

      if ((Teqn(0.00001, tpar) > 0) - (Teqn(0.4, tpar) > 0) == 0) {
        res[i] = uprbd; // no solution within usual timeframe (0.1 hours = 38840psi) so set to upbrd.  otherwise Tc can incorrectly be set to zero.
        tpar[8] = 0.0;
      } else {
        if (T0 < 0.0 || Ts < T0) { // failed during ramp load up to tauc
          res[i] = Ts;
          tpar[8] = 1.0;
        } else {
          if (T0 / Ts < tpar[4]) {  // no damage sustained during constant load
            //res[i] = uprbd;  // this is where we re-do ramp load
            //alpha_T0 = 0;
            alpha_T1 = 0;
            tpar[8] = 0.0;
            //H_T0 = 1;
          } else {  // some damage sustained (or failed) during constant load
            vector<double> C = const_Tc(Ts, T0, tpar);
            //alpha_T0 = C[2];
            alpha_T1 = (C[3] * C[2] - C[0]/C[1] * (  exp(-T1 * C[1]) - C[3] ) ) / exp(-T1 * C[1]);
            tmp = -1 / C[1] * log((C[0] / C[1] * C[3] + C[2] * C[3]) / (1 + C[0] / C[1]));
            tpar[8] = alpha_T1;
            //if (std::isnan(tmp)) res[i] = uprbd;
            //else res[i] = tmp;
          }

          if (alpha_T1 > 1)
            res[i] = tmp;
          else {
            tmp = zero(tpar[4] * Ts + 1e-7, 0.4, t, Teqn2, tpar); // added 1e-7 to prevent numerical error

            res[i] = T1 + tmp;
            if (Teqn2(tmp, tpar) > 0.1) // numerical error no solution
              res[i] = 1e7;

          }

        }
      }
//      // manual checking 
//      cout << setprecision(10) << "a=" << tpar[0] << ", b=" << tpar[1] << ", c=" << tpar[2] << ", n=" << tpar[3] << ", s0=" << tpar[4] << " " ;
//      for (int j = 5; j < 9; j++)
//        cout << tpar[j] << " ";
//      cout << res[i] << endl;
      
    } catch (const exception& e) {
      //			cout << "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
      //			cout << tpar[0] << " " << tpar[1] << " " << tpar[2] << " " << tpar[3] << " " << tpar[4] << " " << tpar[5] << " " << tpar[6] << endl;
    }
  }
  return res;
}

vector<double> summary_stat_v1(vector<double> data, double censored, vector<double> probs) {
  // if Tf > T1, use log(Tf-T1), else use log(Tf)
  // censored is T1  
 
  int n_s = probs.size();
  vector<double> s_obs(n_s, 0);
  
  s_obs = quantile(data, probs);
  
  for (int i = 0; i < n_s; i++)
    s_obs[i] =  (s_obs[i] > censored) ? log(s_obs[i] - censored) : log(s_obs[i]);
  
  return (s_obs);
}


vector<double> summary_stat(vector<double> data, double censored, double k, vector<double> probs) {
  // summary stat for RCR and convert them to strengths.
  
  int n_s = probs.size();
  vector<double> s_obs(n_s, 0);
  
  vector<double> data_mod;
  double count_censored = 0;

  for (int i = 0; i < data.size(); ++i) {
    if (data[i] < censored | std::isnan(data[i])) 
      count_censored += 1;
    else
      data_mod.push_back( (data[i]-censored)*k / 1000.0 ); // strength value
  }
  
  if (data_mod.size() != 0)
    s_obs = quantile(data_mod, probs);
  
  return (s_obs);
}

