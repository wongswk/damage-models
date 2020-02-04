#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <random>
#include <fstream>
#include <sstream>  
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include "model.hpp"
#include "utils.hpp"


using namespace std;


// global variables
const int n_par = 10; // number of parameters
const vector<double> probs = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95}; // quantiles
const int n_s = probs.size(); // number of summary statistics
const vector<double> proposal_scale = {0.01, 0.01, 0.01, 0.01, 0.2, 0.01, 0.01, 0.01, 0.1, 0.01};



// set up random number generators
mt19937 generator(time(0));

int main(int argc, char* argv[])
 {

  if (argc < 6) {
    cout << "Invalid number of samples, burning, thining, delta, and file name(s)" << endl;
    return 0;
  }
  const int m = atoi(argv[1]); // number of ABC samples
  const int burning = atoi(argv[2]);
  const int thining = atoi(argv[3]);
  const double delta = atof(argv[4]);

  int accept = 0;

  uniform_real_distribution<double> unif(0.0, 1.0);
  normal_distribution<double> norm_proposal(0.0, 1);

  ////////////////////////////////////////////////////////////////////////////////

  // read data
  const double k_s = 388440; // MOR PSI per hour
  vector<string> file_name;

  for (int i = 5; i < argc; ++i) {
    file_name.push_back(argv[i]);
  }

  int n_file = file_name.size();

  vector<double> a1(n_file); // MOR PSI
  vector<double> censored(n_file); // 1 year in hours
  vector<double> T0(n_file); // hour
  vector<double> ramprate(n_file);

  string length;
  int pos1, pos2, pos3;
  //vector< vector<double> > load(n_file), data(n_file);
  vector< vector<double> > data(n_file);
  vector<bool> loadType(n_file);  // 0 is 2-phase, 1 is ramp


  for (int i = 0; i < n_file; ++i) {

    string::size_type pos1 = file_name[i].find('_');
    string::size_type pos2 = file_name[i].find('_', pos1 + 1);
    string::size_type pos3 = file_name[i].find('.', pos1 + 1);
    
    //cout << pos1 << " " << pos2 << " " << pos3 << endl;
    
    if (pos2 != -1) { // two-phase
      // extract w0, T0, censored from the file names
      a1[i] = stod(file_name[i].substr(pos1 + 1, pos2 - pos1 - 1));
      T0[i] = a1[i] / k_s;
      length = file_name[i].substr(pos2 + 1, pos3 - pos2 - 1);

      if (length.back() == 'Y') {
        censored[i] = stod(length.substr(0, length.size() - 1)) * 8760;
        if (stod(length.substr(0, length.size() - 1)) == 4)
          censored[i] += 24;  // add a day for leap year
      } else {
        if (length.back() == 'M') {
          censored[i] = stod(length.substr(0, length.size() - 1)) * 8760.0 / 12.0;
        } else {
          cout << "File name format not supported!" << endl;
          return 0;
        }
      }
      loadType[i] = 0;
      ramprate[i] = k_s;
      cout << "input  " << file_name[i] << " 2-phase, tauc = " << a1[i] << ", k_s = " << k_s << ", T1 = " << censored[i] << endl;
    } else { // ramp load (with different rates)
      a1[i] = stod(file_name[i].substr(pos1 + 1, pos3 - pos1 - 1));
      loadType[i] = 1;
      censored[i] = 0.0; //1e9;
      ramprate[i] = a1[i];
      cout << "read " << file_name[i] << " ramp, k = " << a1[i] << endl;
    }
    // read data
    {
      ifstream fin(file_name[i]);
      string line, length;
      int pos;
      if (fin.is_open()) {
        //cout << "Reading file " << file_name[i] << endl;
        while (getline(fin, line)) {

          //string::size_type pos = line.find(',');
          //load[i].push_back(convertToDouble(line.substr(0, pos)));
          data[i].push_back(convertToDouble(line));
        }
        cout << "read " << data[i].size() << " items."  << endl;
      } else {
        cout << "File is not open!" << endl;
      }

    }
  }

  //exit(0);


  ////////////////////////////////////////////////////////////////////////////////


  vector< vector<double> > s_obs(n_file);
  vector<int> N(n_file);
  vector<int> n_c_obs(n_file);

  for (int i = 0; i < n_file; ++i) {
    s_obs[i] = summary_stat(data[i], censored[i], ramprate[i], probs);
    N[i] = data[i].size();
    n_c_obs[i] = count_if(data[i].begin(), data[i].end(), [ = ](double c){return c < censored[i];});
    
    for (int j = 0; j < probs.size(); j++)
      cout << s_obs[i][j] << " ";
    cout << "below T1: " << n_c_obs[i] << endl;
  }

  vector<int> n_obs = N;


  vector< vector<double> > theta(m, vector<double>(n_par));
  vector< vector< vector<double> > > s(m, vector< vector<double> >(n_file));
  vector< vector<double> > newData(n_file);
  vector<double> newtheta(n_par), prevtheta(n_par);
  vector< vector<double> > s_simu(n_file, vector<double>(n_s)), prev_s(n_file, vector<double>(n_s));
  double u;
  double log_alpha;
  double log_den = 0;
  double log_num = 0;
  vector<int> n_censored(n_file);



  // generate random samples
  theta[0] = {-8.33547, 0.527351, 3.77781, 0.0567059, -24.2165, 0.0897585, -0.847272, 0.405729, -0.737646, 0.179287}; // the initial value of theta

  log_den = log_prior(theta[0], n_par);
  for (int f = 0; f < n_file; ++f) {
    if (loadType[f])
      newData[f] = CANR(theta[0], N[f], a1[f], k_s);
    else
      newData[f] = CANL(theta[0], N[f], a1[f], k_s, censored[f]);
    s[0][f] = summary_stat(newData[f], censored[f], ramprate[f], probs);
    n_censored[f] = count_if(newData[f].begin(), newData[f].end(), [ = ](double c){return c < censored[f];});
    //log of the denominator in alpha
    log_den += dnorm(dist(s[0][f], s_obs[f]) / delta, 0.0, 1.0, true);
    if (n_c_obs[f] > 0) {
      log_den += n_c_obs[f] * log((double) n_censored[f] / (double) N[f]) +
              (n_obs[f] - n_c_obs[f]) * log(1 - (double) n_censored[f] / (double) N[f]);
    }
  }



  // burning	
  prevtheta = theta[0];
  prev_s = s[0];
  for (int i = 1; i < burning; ++i) {

    // sample parameters from the proposal distribution
    for (int l = 0; l < n_par; ++l) {
      newtheta[l] = prevtheta[l] + proposal_scale[l] * norm_proposal(generator);
    }


    if (newtheta[1] < 0 || newtheta[3] < 0 || newtheta[5] < 0 || newtheta[7] < 0 || newtheta[9] < 0) {
      continue;
    }

    // generate data
    log_num = log_prior(newtheta, n_par);
    for (int f = 0; f < n_file; ++f) {
      if (loadType[f])
        newData[f] = CANR(newtheta, N[f], a1[f], k_s);
      else
        newData[f] = CANL(newtheta, N[f], a1[f], k_s, censored[f]);
      
      s_simu[f] = summary_stat(newData[f], censored[f], ramprate[f], probs);
      
//      for (int j = 0; j < probs.size(); j++)
//        cout << s_simu[f][j] << " ";
//      cout << endl;
      n_censored[f] = count_if(newData[f].begin(), newData[f].end(), [ = ](double c){return c < censored[f];});
      //log of the numerator in alpha
      log_num += dnorm(dist(s_simu[f], s_obs[f]) / delta, 0.0, 1.0, true);

      if (n_c_obs[f] > 0) {
        log_num += n_c_obs[f] * log((double) n_censored[f] / (double) N[f]) +
                (n_obs[f] - n_c_obs[f]) * log(1 - (double) n_censored[f] / (double) N[f]);
      }
      if (n_c_obs[f] == 0 & n_censored[f] > 0) {
        log_num -= 10000.0;
      }
    }

    // accept/reject step
    u = unif(generator);

    log_alpha = log_num - log_den;
    if (log(u) < log_alpha &
            !std::isnan(log_alpha) &
            all_of(newData.begin(), newData.end(), [](vector<double> d) {
              return all_of(d.begin(), d.end(), [](double dd) {
                return dd > 0;
              });
            })) {

    prevtheta = newtheta;
    log_den = log_num;
    prev_s = s_simu;
  }
  }


  theta[0] = prevtheta;
  s[0] = prev_s;

  for (int i = 1; i < m; ++i) {

    for (int j = 0; j < thining; ++j) {

      // sample parameters from the proposal distribution
      for (int l = 0; l < n_par; ++l) {
        newtheta[l] = prevtheta[l] + proposal_scale[l] * norm_proposal(generator);
      }

      if (newtheta[1] < 0 || newtheta[3] < 0 || newtheta[5] < 0 || newtheta[7] < 0 || newtheta[9] < 0) {
        continue;
      }


      // generate data
      log_num = log_prior(newtheta, n_par);
      for (int f = 0; f < n_file; ++f) {
        if (loadType[f])
          newData[f] = CANR(newtheta, N[f], a1[f], k_s);
        else
          newData[f] = CANL(newtheta, N[f], a1[f], k_s, censored[f]);
        s_simu[f] = summary_stat(newData[f], censored[f], ramprate[f], probs);
        n_censored[f] = count_if(newData[f].begin(), newData[f].end(), [ = ](double c){return c < censored[f];});
        //log of the numerator in alpha
        log_num += dnorm(dist(s_simu[f], s_obs[f]) / delta, 0.0, 1.0, true);
        cout << dist(s_simu[f], s_obs[f]) << "(" << n_censored[f] << ") ";
        if (n_c_obs[f] > 0) {
          log_num += n_c_obs[f] * log((double) n_censored[f] / (double) N[f]) +
                  (n_obs[f] - n_c_obs[f]) * log(1 - (double) n_censored[f] / (double) N[f]);
        }
        if (n_c_obs[f] == 0 & n_censored[f] > 0) {
          log_num -= 10000.0;
        }

      }
      cout << log_num << " ";
      
      // accept/reject step
      u = unif(generator);

      log_alpha = log_num - log_den;
      if (log(u) < log_alpha &
              !std::isnan(log_alpha) &
              all_of(newData.begin(), newData.end(), [](vector<double> d) {
                return all_of(d.begin(), d.end(), [](double dd) {
                  return dd > 0;
                });
              })) {

      prevtheta = newtheta;
      log_den = log_num;
      prev_s = s_simu;
    }

      if (log(u) < log_alpha &
              !std::isnan(log_alpha) &
              all_of(newData.begin(), newData.end(), [](vector<double> d) {
                return all_of(d.begin(), d.end(), [](double dd) {
                  return dd > 0;
                });
              })) {

      prevtheta = newtheta;
      log_den = log_num;
      prev_s = s_simu;
      cout << "ACCEPT ";
      accept++;
    }
              cout << endl;
    }
    theta[i] = prevtheta;
    s[i] = prev_s;
  }

  // the output file names depend on delta
  stringstream theta_name;
  //stringstream stat_name;
  stringstream accept_name;

  theta_name << "theta_" << delta << ".csv";
  //stat_name << "stat_" << delta << ".csv";
  accept_name << "accept_" << delta << ".csv";


  // output simulation result  
  ofstream fout(theta_name.str());

  if (fout.is_open()) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n_par - 1; j++) {
        fout << theta[i][j] << ",";
      }
      fout << theta[i][n_par - 1] << endl;
    }
  } else {
    cout << "File is not open!" << endl;
  }

  ofstream fout3(accept_name.str());

  if (fout3.is_open()) {
    fout3 << (double) accept / (double) (m * thining) << endl;
  } else {
    cout << "File is not open!" << endl;
  }

  return 0;
}

