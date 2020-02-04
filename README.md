# damage-models
Calibrating wood products for load duration and rate using US, Canadian, and gamma process models

## Overview

This repository contains programs to fit the US, Canadian, and gamma process models and perform reliability analysis.  Code is composed of R scripts and C++ source files.

### Download and compile the C++ programs for Linux:
```
git clone https://github.com/wongswk/damage-models.git
cd damage-models/src
make all
```
### Example data

Two example datasets are included in "data/" for usage demonstrations in each model below.

1. ramp_388440.txt:  Ramp load with rate k = 388440 psi/hour

2. RCR_4500_1Y.txt:  Ramp-constant-ramp load, i.e., constant load 4500psi for duration 1 year, followed by a ramp load test to failure of survivors at the standard rate

## US (Gerhard's) damage model

### Fitting the model

First, edit 'R/loadData.R' as needed to include the datasets and settings for the experiment.  The included file is set up to handle the example datasets above, and contains instructions for how to make modifications.

In R, set the working directory to the base folder of this repository (damage-models)

Example:
```
source("R/fit_US.R")
```
This fits the US model and writes out the parameter estimates to 'US_theta.csv'.  The median short-term strength can be set in Line 3 of 'R/fit_US.R'.

### Reliability analysis

Reliability analysis is performed by the compiled C++ program 'USADM'

Example:
```
./USADM US_theta.csv 1.2  

```
Using the values of theta listed in "US_theta.csv", this program will estimate the probability of failure after 50 years based on a large
	number of replications with the stochastic load profile using _phi=1.2_.

General usage:

USADM **theta_file_name** **phi** 

- The parameter vectors to use for reliability analysis are supplied via **theta_file_name**, usually the output of the parameter estimation done in model fitting above.
	The thetas will then be used to solve the US ADM for time-to-failure T_f and
	probability of failure P_f. T_f and P_f for each theta will be written into files with suffix “csv”.
- Output files: probability of failure with and without DOL effect.

There are some additional parameters in USADM (solveUSODE.cpp) that can be modified:

- n_per_theta:    number of load profile replications to use for estimation of P_f.  Default is 100000.
- t_end:        reference timeframe (in hours).  Default is 50 years.
- dt:        time-step used for numerical integration.  Default is 100 hours.
- The current load profile is the residential load profile. Other load profiles can be added by creating a new class.


## Canadian (Foschi's) damage model

### Fitting the model

Model fitting via ABC-MCMC is performed by the compiled C++ program 'DOLmod_MCMC'

Example:
```
./DOLmod_MCMC 50 1000 100 2.0 data/ramp_388440.txt data/RCR_4500_1Y.txt
```
This generates 50 MCMC samples (thinned by every 100 iterations), after a burn-in of 1000 iterations, using an ABC delta tolerance of 2.0.  The data fitted are the files data/ramp_388440.txt and data/RCR_4500_1Y.txt

General usage:

DOLmod_MCMC **m** **burnin** **thining** **delta** **data**

- The **data** can consist of multiple file names, each named following the format in the example datasets
- Output files: theta.csv containing _m_ theta vectors and accept.csv containing the acceptance rate.


### Reliability analysis

Reliability analysis is performed by the compiled C++ program 'CanADM'

Example:
```
./CanADM theta_2.csv 1.2  

```
For each sample of theta in "theta_2.csv", this program will estimate the probability of failure after 50 years based on a large
	number of replications with the stochastic load profile using _phi=1.2_.

General usage:

CanADM **theta_file_name** **phi** 

- The parameter vectors to use for reliability analysis are supplied via **theta_file_name**, usually the output of 'DOLmod_MCMC' program above.
	The thetas will then be used to solve the Canadian ADM for time-to-failure T_f and
	probability of failure P_f. T_f and P_f for each theta will be written into files with suffix “csv”.
- Output files: probability of failure with and without DOL effect.

There are some additional parameters in CanADM (solveCanODE.cpp) that can be modified:

- n_per_theta:    number of load profile replications to use for estimation of P_f.  Default is 100000.
- t_end:        reference timeframe (in hours).  Default is 50 years.
- dt:        time-step used for numerical integration.  Default is 100 hours.
- The current load profile is the residential load profile. Other load profiles can be added by creating a new class.


## Gamma process model

First, edit 'R/loadData.R' as needed to include the datasets and settings for the experiment.  The included file is set up to handle the example datasets above, and contains instructions for how to make modifications.  'tau_max' should be set to the maximum strength of the lumber population.

In R, set the working directory to the base folder of this repository (damage-models)

### Fitting the model

Example:
```
source("R/fit_GP.R")
```
This fits the gamma process model and writes out the sampled parameters to 'fit-gp/1.rda'.

The main settings in the R script (R/fit_GP.R) that can be modified:

- The gamma process model file, default is **R/gamma-2breaks.R**, which allows for two breakpoints in the power law.  Model files for 0, 1, 2, and 3 breakpoints are included in the package (R/gamma-0breaks.R, R/gamma-1breaks.R, R/gamma-2breaks.R, R/gamma-3breaks.R)
- n.iter:  number of MCMC iterations to run
- nprocs: number of parallel threads (CPU cores) to use
- Other settings for parallel tempering that can be modified are described within the R script

### Reliability analysis

Begin by generating a large number of stochastic load profiles.

Example:
```
source("R/setup-reliability.R")
```
The number load profile replications to generate (default: 100000), and parameters associated with the load profiles can be modified within that R script.

Now calculate probability of failure for each profile, using parameters obtained from model fitting.

Example:
```
source("R/calc-reliability-GP.R")
```
The results are saved as a matrix of failure probabilities **failmat**.

The default model file is **R/gamma-2breaks.R**, and should match the model used for fitting.  The performance factor **phi** and reference time frame (default 50 years) can be set within the R script.

## References

Some code in this repository is adapted from the methods described in these papers:

Samuel WK Wong and James V Zidek (2019). The duration of load effect in lumber as stochastic degradation. IEEE Transactions on Reliability, 68(2), 410-419.

Chun-Hao Yang, James V Zidek, Samuel WK Wong (2019). Bayesian analysis of accumulated damage models in lumber reliability. Technometrics, 61(2), 233-245.

Questions?  Contact: samuel.wong@uwaterloo.ca

