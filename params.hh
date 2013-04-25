//file "params.hh"

// membrane time constant tau_m=10 msec = 0.01 sec

#ifdef DEFINE_GLOBAL
#	define EXTERN 
#else
#	define EXTERN extern
#endif

EXTERN unsigned int C_ndt, C_ndt_mu, C_N_trials, C_N_neuron, C_N_Gen, C_rho_k_max, C_tau_r__dt, C_size_powspe;
EXTERN double C_rate, C_dt, C_eps, C_tau_r, C_tau_s, C_tau_ou, C_D_ou, C_T_max, C_T_mu, C_tol_mu, C_df, M2_PI;

#undef EXTERN
/*
double C_rate= 1.; // firing rate [1/tau_m]=10^{2} [Hz]
const double C_dt= .001; // timestep [tau_m]=10^{-2}[sec]
const unsigned int C_ndt=1048576;//8388608;//4194304;//2097152;//1048576;//524288;//262144;//131072;//65536; // number of timesteps in simulation
const unsigned int C_ndt_mu=C_ndt/2;//32768;//16384; // number of timesteps in bisection
const unsigned int C_N_trials=10; // number of trials per generation
double C_eps=1.; // weight of Input-Current
const double C_tau_r=.1; // absolute refractory period [tau_m]=10^{-2}[sec]
const unsigned int C_N_neuron=100 ; // number of neurons per generation
const unsigned int C_N_Gen=1 ; // number of layers in the feedforward-network (generations)

const double C_tau_s=1. ; // time constant for synaptic filtering - 2-4 msec (cut-off frequency for powerspectrum at 1/tau_s) [tau_m]=10^{-2}[sec]

const unsigned int C_rho_k_max=11; // max. lag number for serial correlation coefficient

const double C_tau=.5; // time constant for ornstein-uhlenbeck-powerspectrum
double C_D=C_rate/C_tau/4; // diffusion coefficient for ornstein-uhlenbeck-powerspectrum

const unsigned int C_tau_r__dt=C_tau_r/C_dt; // absolute refractory period [#timesteps]
const double C_T_max=C_dt*(C_ndt-1); // simulation-time [tau_m]=10^{-2}[sec]
const double C_T_mu=C_dt*(C_ndt_mu-1); // bisection simulation-time [tau_m]=10^{-2}[sec]
const double C_tol_mu=1./C_T_mu; // tolerance for bisection
const unsigned int C_size_powspe=C_ndt/2+1; // size of vector for powerspectra

const double M2_PI=2*M_PI; // 2x pi
const double MSQR_PI=pow(M_PI,2); // pi^2

/* some values to know:

	f_max = 1/(2*dt) = 500 * 100 Hz
	df = 1/(N*dt) = 0.00023 * 100 Hz for simulation
	df = 1/(N_mu*dt) = 0.0019 * 100 Hz for mutest

*/
