// file "main.cc"

#define DEFINE_GLOBAL
#include "params.hh"
#undef DEFINE_GLOBAL
#include "functs.hh"
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <omp.h>
	using namespace std;

int main(int argc, char *argv[])
{
// initialization of global variables
	C_rate= 1.; // firing rate [1/tau_m]=10^{2} [Hz]
	C_dt= .001; // timestep [tau_m]=10^{-2}[sec]
	C_ndt=1048576;//8388608;//4194304;//2097152;//1048576;//524288;//262144;//131072;//65536; // number of timesteps in simulation
	C_ndt_mu=C_ndt/2;//32768;//16384; // number of timesteps in bisection
	C_N_trials=10; // number of trials per generation
	C_eps=1.; // weight of Input-Current
	C_tau_r=.1; // absolute refractory period [tau_m]=10^{-2}[sec]
	C_N_neuron=100 ; // number of neurons per generation
	C_N_Gen=2 ; // number of layers in the feedforward-network (generations)

	C_tau_s=1. ; // time constant for synaptic filtering - 2-4 msec (cut-off frequency for powerspectrum at 1/tau_s) [tau_m]=10^{-2}[sec]

	C_rho_k_max=11; // max. lag number for serial correlation coefficient

	C_tau_ou=.5; // time constant for ornstein-uhlenbeck-powerspectrum
	C_D_ou=C_rate/C_tau_ou/4.; // diffusion coefficient for ornstein-uhlenbeck-powerspectrum

	C_tau_r__dt=C_tau_r/C_dt; // absolute refractory period [#timesteps]
	C_T_max=C_dt*(C_ndt-1.); // simulation-time [tau_m]=10^{-2}[sec]
	C_T_mu=C_dt*(C_ndt_mu-1.); // bisection simulation-time [tau_m]=10^{-2}[sec]
	C_tol_mu=1./C_T_mu; // tolerance for bisection
	C_size_powspe=C_ndt/2+1; // size of vector for powerspectra
	C_df=1./(C_ndt*C_dt);

	M2_PI=2*M_PI; // 2x pi


// variables and constants for calculation
	double *powspe_old, *powspe_new, *S_temp, *I_input;
	double mu_gen[C_N_trials], rate_gen[C_N_trials], sigma_gen[C_N_trials], CV_gen[C_N_trials], rho_k_gen[C_rho_k_max*C_N_trials], mu;
	double mu_avg, rate_avg, sigma_avg, CV_avg, rho_avg[C_rho_k_max], rho_var[C_rho_k_max];
	
	dzeros(rho_avg, C_rho_k_max); dzeros(rho_var, C_rho_k_max);

// variables for initialization of random number generators	
	const time_t now=time(NULL);

// variables for saving data
	stringstream buffer;
	char date[18];
	strftime(date,18, "%Y-%m-%d_%H-%M",localtime(&now));
	stringstream filename;
	string filename_tmp;

//	fftw_init_threads();
//	fftw_plan_with_nthreads(omp_get_max_threads());

// create plan for fft in "calc_I_diff"
	I_input=new double[C_ndt];
	fftw_plan plan_I_input = fftw_plan_r2r_1d(C_ndt, I_input, I_input, FFTW_HC2R, FFTW_MEASURE | FFTW_DESTROY_INPUT);
	char* fft_wisdom_I_input=fftw_export_wisdom_to_string();
	delete[] I_input;
	fftw_destroy_plan(plan_I_input);

// generate start-powerspectrum and save it
	powspe_old=new double[C_size_powspe];
	ornstein_PS(powspe_old); 
	safe_powspe(powspe_old, rho_avg, rho_var, C_rate, C_tau_ou, C_D_ou, log(-1), log(-1), log(-1), log(-1), log(-1), 0, date);
//	whitenoise_PS(powspe_old,C_rate);
//	safe_powspe(powspe_old, rho_avg, rho_var, C_rate, log(-1), log(-1), log(-1), log(-1), log(-1), log(-1), log(-1), 0, date);
		
	powspe_new=new double[C_size_powspe];
	dzeros(powspe_new, C_size_powspe);

// ----------------------------------------------------------------- //
// Startpoint: Simulation
	for (unsigned int i_gen = 1; i_gen <= C_N_Gen; i_gen++){
/* T */		cout << "Gen: " << i_gen << endl;

		for (unsigned int i_trial = 0; i_trial < C_N_trials; i_trial++){
/* T */			cout << "Neuron: " << i_trial << endl;

			ISI interval(C_T_max);
			I_input=new double[C_ndt];
	
		// generate Diffusion-Current
			fftw_import_wisdom_from_string(fft_wisdom_I_input);
			plan_I_input = fftw_plan_r2r_1d(C_ndt, I_input, I_input, FFTW_HC2R, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
			calc_I_diff(I_input, powspe_old, plan_I_input, ((i_gen-1)*C_N_trials+i_trial+static_cast<unsigned long>(now)));
			fftw_destroy_plan(plan_I_input);

		// define mu (and eps_diff) + create ISI-train 
				mu = mutest(I_input);
				interval.lif_neuron(mu, I_input, C_ndt);

			delete[] I_input;

		// calculate Powerspectrum - normalized to the rate 
			S_temp= new double[C_size_powspe];
			powerspectrum(S_temp, interval);

			for (unsigned int i_S = 0; i_S < C_size_powspe; i_S++){
				powspe_new[i_S]+=S_temp[i_S]/C_N_trials;
			}			
			delete[] S_temp;

		// calculate averages of rate and mu 
			rate_gen[i_trial]=interval.rate();
			sigma_gen[i_trial]=interval.var();
			CV_gen[i_trial]=interval.CV();
			mu_gen[i_trial]= mu;
			interval.calc_rho_k(rho_k_gen, i_trial, interval);

		}

	// calculating mean and variance of serial correlation coefficient rho_k
		for (unsigned int k = 0; k < C_rho_k_max; k++){
			rho_avg[k]=0;	
			rho_var[k]=0;
			for (unsigned int i = 0; i < C_N_trials; i++){
				rho_avg[k]+=rho_k_gen[k+i*C_rho_k_max]/C_N_trials;
				rho_var[k]+=(rho_avg[k]-rho_k_gen[k+i*C_rho_k_max])*(rho_avg[k]-rho_k_gen[k+i*C_rho_k_max])/(C_N_trials-1);
			}
		}

	// saving powerspectrum, rate, sigma, mu and constants to file
		rate_avg=average(rate_gen,C_N_trials);
		sigma_avg=average(sigma_gen,C_N_trials);
		CV_avg=average(CV_gen,C_N_trials);
		mu_avg=average(mu_gen,C_N_trials);

		safe_powspe(powspe_new, rho_avg, rho_var, rate_avg, sigma_avg, CV_avg, mu_avg, variance(rate_gen,C_N_trials,rate_avg), variance(sigma_gen,C_N_trials,sigma_avg), variance(CV_gen,C_N_trials,CV_avg), variance(mu_gen,C_N_trials,mu_avg), i_gen, date);
	
		cpNcl_d_arr(powspe_new, powspe_old, C_size_powspe);
	}

// free all
	delete[] powspe_new;
	delete[] powspe_old;
	fftw_forget_wisdom();
	fftw_cleanup();//_threads();

	return 0;
}
