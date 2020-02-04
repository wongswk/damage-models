#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <stdlib.h>
#include <math.h>
#include <iostream> // for cin and cout
#include <fstream>
#include <sstream>  
#include <vector> // for vector
#include <map>
#include <string>
#include <omp.h>
#include <random> // for distributions
#include "boost/math/special_functions/gamma.hpp" // for the incomplete gamma function
#include "boost/numeric/odeint.hpp" // for the ode solver
#include "boost/math/tools/roots.hpp" // for the root finding

using namespace std;
using namespace boost::numeric::odeint;
using boost::math::tgamma_lower;
using namespace boost::math::tools; // for bisect

////////////////////////////////////////////////////////////////////////////////

// setup parameters
const double mu = 1;
const double k = 388440;
const double T0 = 4500/k;
const int n_per_theta = 100000;
const eps_tolerance<double> tol(16);
const double t_start = 0;
const double t_end = 438300; // hours in 50 years
double dt = 100;  // integation step size

map<string, double> load_profile_par{
	{"phi", 1},
	{"load_s_shape", 3.122}, 
	{"load_s_scale", 0.0481}, 
	{"load_p_shape", 0.826}, 
	{"load_p_scale", 0.1023}, 
	{"R0", 3000},
	{"alpha_d", 1.25},
	{"alpha_l", 1.5},
	{"load_d_mean", 1.05},
	{"load_d_sd", 0.1},
	{"gamma", 0.25},
	{"mean_Ts", 10},
	{"mean_Te", 1},
	{"mean_Tp", 0.03835},
	{"N", 100} };

// set up random number generators
mt19937 generator(0);

////////////////////////////////////////////////////////////////////////////////

double stepfun(double x, vector<double> x_vec, vector<double> y_vec){
	if(y_vec.size() - x_vec.size() != 1){
		cout << "The length of y_vec must be greater than the length of x_vec by 1." << endl;
		return(0);
	}
	int n = x_vec.size();
	if(x < x_vec[0]) return y_vec[0];
	else if(x >= x_vec[n-1]) return y_vec[n];
	else{
		for(int i = 1; i < n; ++i){
			if(x < x_vec[i] && x >= x_vec[i-1]) return y_vec[i];
		}
	}
}

class load_profile {
	
	double load_p_shape, load_p_scale, load_s_shape, load_s_scale, phi, R0, alpha_d, alpha_l,
		mean_Ts, mean_Te, mean_Tp, D_d, load_d_mean, load_d_sd, gamma;
	int N = 250;


public:
	vector<double> T_s = vector<double>(N);
	vector<double> T_e = vector<double>(N);
	vector<double> load_s = vector<double>(N+1);
	vector<double> load_e = vector<double>(N+1);

	load_profile(map<string, double> load_profile_par):
		load_s_shape(load_profile_par["load_s_shape"]),
		load_s_scale(load_profile_par["load_s_scale"]),
		load_p_shape(load_profile_par["load_p_shape"]),
		load_p_scale(load_profile_par["load_p_scale"]),
		phi(load_profile_par["phi"]),
		alpha_d(load_profile_par["alpha_d"]),
		alpha_l(load_profile_par["alpha_l"]),
		R0(load_profile_par["R0"]),
		load_d_mean(load_profile_par["load_d_mean"]),
		load_d_sd(load_profile_par["load_d_sd"]),
		gamma(load_profile_par["gamma"]),
		mean_Ts(load_profile_par["mean_Ts"]),
		mean_Te(load_profile_par["mean_Te"]),
		mean_Tp(load_profile_par["mean_Tp"]) 
		{
			normal_distribution<double> norm_d(load_d_mean, load_d_sd);
			gamma_distribution<double> gamma_s(load_s_shape, load_s_scale);
			gamma_distribution<double> gamma_p(load_p_shape, load_p_scale);
			exponential_distribution<double> exp_Ts(1.0/mean_Ts);
			exponential_distribution<double> exp_Te(1.0/mean_Te);
			exponential_distribution<double> exp_Tp(1.0/mean_Tp);
			D_d = norm_d(generator);	
			T_s[0] = exp_Ts(generator);			
			T_e[0] = exp_Te(generator);
			load_s[0] = gamma_s(generator);
			load_e[0] = 0;
			for(int i = 1; i < N; ++i){
				T_s[i] = T_s[i-1] + exp_Ts(generator);
				load_s[i] = gamma_s(generator);
				if(i % 2 != 0){
					T_e[i] = T_e[i-1] + exp_Tp(generator);
					load_e[i] = gamma_p(generator);
				}
				else{
					T_e[i] = T_e[i-1] + exp_Te(generator);
					load_e[i] = 0;				
				}
			}
			load_s[N] = gamma_s(generator);
			if(load_e[N-1] == 0) load_e[N] = 0;
			else load_e[N] = gamma_p(generator);
		}
 

	double operator() (double t){
		t = t / 8760; // express t in years
		double D_s = stepfun(t, T_s, load_s);
		double D_e = stepfun(t, T_e, load_e);
		double load_t = phi * R0 * (gamma * D_d + D_e + D_s)/(gamma * alpha_d + alpha_l);
		return load_t;
	}

	vector<double> max(){
		double max_load;
		vector<double> res(2,0);
		for(double t = t_start; t < t_end; t += dt){
			max_load = (*this)(t);
			if(max_load > res[1]){
				res[0] = t;
				res[1] = max_load;			
			}	
		}
		return(res);
	}
	
	void write(ofstream& file){
		if(file.is_open()){
			file << "T_s, load_s, T_e, load_e" << endl;
			for(int i = 0; i < T_s.size(); ++i){
				file << T_s[i] << "," << load_s[i] << "," << T_e[i] << "," << load_e[i] << endl;
			}	
			file << " ," << load_s[T_s.size()] << ", ," << load_e[T_s.size()] << endl;
		}
		else{
			cout << "File is not open!" << endl;
		}
	}

};

class constant_load{

	double load_level;	

public:
	constant_load(double l): load_level(l) {};

	double operator() (double t){
		// t is in hour
		double T0 = load_level/k;
		if(t <= T0) return k*t;
		else return load_level;
	}

	vector<double> max(){
		double T0 = load_level/k;
		vector<double> res(2,0);
		res[0] = T0;
		res[1] = load_level;
		return(res);
	}

};

template <class T>
void CanADM(const double &x , double &dxdt , const double t, vector<double> parameter, double tau_s, T tau){
	
	double a, b ,c ,n, s0, mu;

	a = parameter[0];
	b = parameter[1];
	c = parameter[2];
	n = parameter[3]; 
	s0 = parameter[4]; 
	mu = parameter[6]; 

	dxdt = 1/mu * ( pow(a * tau_s * max((tau(t)/tau_s - s0), 0.0), b) + ( pow(c * tau_s * max((tau(t)/tau_s - s0), 0.0), n) * x));
}

double Teqn(double Ts, vector<double> params) {
	double A = params[0];
	double b = params[1];
	double C = params[2];
	double n = params[3];
	double s0 = params[4];
	double k = params[5];
	double mu = params[6];
  
	double As = A*k*Ts;
	double Cs = C*k*Ts;
  
	double intfac = 1/mu * pow(Cs, n)* Ts/(n+1) * pow(1-s0,n+1);
	double bnfr = (b+1)/(n+1);
  
	double result = pow(As,b) / pow(Cs, n*bnfr) * pow(mu/Ts * (n+1),(b-n)/(n+1)) * tgamma_lower(bnfr, intfac) - exp(-intfac);
      
	return result;
}


class BadConversion : public std::runtime_error {
public:
  BadConversion(const std::string& s)
    : std::runtime_error(s)
    { }
};

double convertToDouble(const string& s)
{
    istringstream i(s);
    double x;
    char c;
    if (!(i >> x))
        throw BadConversion("convertToDouble(\"" + s + "\")");
    return x;
}

vector< vector<double> > readCSV_to_vec(ifstream& file){
	vector< vector<double> > result;
	vector<double> row;
	double p;
	string line;
    string cell;

	while(getline(file, line)){
		
		stringstream lineStream(line);

		while(getline(lineStream, cell, ','))
		{
			p = convertToDouble(cell);
		    row.push_back(p);
		}
		result.push_back(row);
		row.clear();
	}
	return result;
}

vector<double> generate(vector<double> &theta){
	// Generate random variates
	normal_distribution<double> normA(theta[0], theta[1]);
	normal_distribution<double> normb(theta[2], theta[3]);
	normal_distribution<double> normC(theta[4], theta[5]);
	normal_distribution<double> normn(theta[6], theta[7]);
	normal_distribution<double> norms0(theta[8], theta[9]);

	vector<double> pars(7);
	pars[0] = exp(normA(generator));
	pars[1] = exp(normb(generator));
	pars[2] = exp(normC(generator));
	pars[3] = exp(normn(generator));
	double tempx = exp(norms0(generator));
	pars[4] = tempx / (1.0+tempx);  
	pars[5] = k;
	pars[6] = mu;

	return pars;
}



////////////////////////////////////////////////////////////////////////////////


// main function
int main(int argc, char* argv[]){

	double Ts, tau_s, x, t;
	double x_max = 1;


	if(argc == 3) load_profile_par["phi"] = atof(argv[2]);
	cout << "phi is set to be " << load_profile_par["phi"] << endl;
        
        //cout << tol << endl;  

	// read theta from file
	ifstream fin(argv[1]);
	if(fin.is_open()){
		vector< vector<double> > theta = readCSV_to_vec(fin);
		int ntheta = theta.size();

		vector<int> nFail_per_theta(ntheta, 0);
		vector<int> nFail_per_theta_noDOL(ntheta, 0);
		vector<double> prop_Fail(ntheta, 0);
		vector<double> prop_Fail_noDOL(ntheta, 0);
		vector<double> time_to_failure(ntheta*n_per_theta);	

		vector<double> pars;	

		cout << "Read " << ntheta << " thetas!" << endl;
		cout << "Start simulation!" << endl;

		

		omp_set_num_threads(16);
		#pragma omp parallel for schedule(dynamic) private(t, x, pars, Ts, tau_s)		
		for(int i = 0; i < ntheta; ++i){
			int tmp = 0;
			int tmp_noDOL = 0;

			
			adams_bashforth_moulton< 5 , double > stepper;

			for(int j = 0; j < n_per_theta; ++j){

				
				
				pars = generate(theta[i]);

				load_profile tau(load_profile_par);
				//constant_load tau(4500);

				try{
					pair<double, double> r = bisect([pars](double t) {return Teqn(t, pars);}, 0.001, 0.5, tol);
					Ts = r.first + (r.second - r.first)/2;
                                  
				}
				catch(const boost::math::evaluation_error& e){
					time_to_failure[i * n_per_theta + j] = t_end/8760.0;
					continue;
				}

				tau_s = Ts * k;

				vector<double> max_t = tau.max();

				if(max_t[1] < tau_s * pars[4]){
					time_to_failure[i * n_per_theta + j] = t_end/8760.0;
				}
				else{
					
					// no DOL
					if(max_t[1] > tau_s){
    					tmp_noDOL += 1;
						tmp += 1;
						continue;
					} 
					
					// solve the ODE 		
					//CanADM model(pars, tau_s, tau);
					auto model = [&pars, &tau_s, &tau](const double &x , double &dxdt , const double t){
						CanADM(x , dxdt , t, pars, tau_s, tau);	
					};			


					x = 0; // initial state
					t = t_start;
					//times.push_back(t);
					//x_vec.push_back(x);
					stepper.initialize(model, x , t , dt );
					t += dt;
					while(t < t_end){
						stepper.do_step(model, x, t, dt);
						//times.push_back(t);
						//x_vec.push_back(x);
						if( x > x_max || std::isnan(x)) // if x > 1, stop calculation
						{
							tmp += 1;
							break;
						}
						t += dt;
					}

					time_to_failure[i * n_per_theta + j] = t/8760.0;


				}
				
				// the ode can be solved simply with the following; but we want it stops early when some criterion is met
				// so we have to solve it manually.

    			//integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) ,
                //    model , x , t_start , t_end , dt, push_back_state_and_time(x_vec, times));

			}
			nFail_per_theta[i] = tmp;
			nFail_per_theta_noDOL[i] = tmp_noDOL;
			prop_Fail[i] = (double)nFail_per_theta[i] / (double) n_per_theta;
			prop_Fail_noDOL[i] = (double)nFail_per_theta_noDOL[i] / (double) n_per_theta;

			cout << i << " " << tmp << " " << tmp_noDOL << endl;

		}

		// the output file names depend on phi
		stringstream prob_name;
		stringstream prob_noDOL_name;
		stringstream time_name;

		prob_name << "prob_" << load_profile_par["phi"] << ".csv";
		prob_noDOL_name << "prob_noDOL_" << load_profile_par["phi"] << ".csv";	
		time_name << "time_" << load_profile_par["phi"] << ".csv";	

		
		// output time_to_failure
		ofstream file2(time_name.str());
		if(file2.is_open()){
			for(int i = 0; i < time_to_failure.size(); ++i){
				file2 << time_to_failure[i] << endl;
			}
		}
		else{
			cout << "File is not open!" << endl;
		}
		
		
		ofstream file3(prob_name.str());
		if(file3.is_open()){
			for(int i = 0; i < prop_Fail.size(); ++i){
				file3 << prop_Fail[i] << endl;
			}
		}
		else{
			cout << "File is not open!" << endl;
		}
		

		ofstream file4(prob_noDOL_name.str());
		if(file4.is_open()){
			for(int i = 0; i < prop_Fail.size(); ++i){
				file4 << prop_Fail_noDOL[i] << endl;
			}
		}
		else{
			cout << "File is not open!" << endl;
		}

	
	}
	else{
		cout << "Please specify the file of theta!" << endl;
	}

	return(0);
}

////////////////////////////////////////////////////////////////////////////////



