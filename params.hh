//file "params.hh"

// membrane time constant tau_m=10 msec = 0.01 sec

const double C_rate= .1; // firing rate [1/tau_m]=10^{2} [Hz]
const double C_cv=.3 ; // coefficient of variation
const double C_dt= .001; // timestep [tau_m]=10^{-2}[sec]
const unsigned int C_ndt=1048576;//524288;//262144;//131072;//65536; // number of timesteps in simulation
const unsigned int C_ndt_mu=C_ndt/2;//32768;//16384; // number of timesteps in bisection
const double C_eps=.3; // weight of Input-Current
const double C_tau_r=.1; // absolute refractory period [tau_m]=10^{-2}[sec]
const unsigned int C_N_neuron=5 ; // number of neurons per generation
const unsigned int C_N_Gen=1 ; // number of layers in the feedforward-network (generations)

const double C_tau_s=.3 ; // time constant for synaptic filtering - 2-4 msec (cut-off frequency for powerspectrum at 1/tau_s) [tau_m]=10^{-2}[sec]

const unsigned int C_rho_k_max=11; // max. lag number for serial correlation coefficient

const double C_D=5.; // diffusion coefficient for ornstein-uhlenbeck-powerspectrum
const double C_tau=.5; // time constant for ornstein-uhlenbeck-powerspectrum

const unsigned int C_tau_r__dt=C_tau_r/C_dt; // absolute refractory period [#timesteps]
const double C_T_max=C_dt*(C_ndt-1); // simulation-time [tau_m]=10^{-2}[sec]
const double C_T_mu=C_dt*(C_ndt_mu-1); // bisection simulation-time [tau_m]=10^{-2}[sec]
const double C_tol_mu=1./C_T_mu; // tolerance for bisection
const unsigned int C_size_powspe=C_ndt/2+1; // size of vector for powerspectra

const double M2_PI=2*M_PI; // 2x pi
const double MSQR_PI=pow(M_PI,2); // pi^2

/* some values to know:

	f_max = 1/(2*dt) = 500 * 100 Hz
	df = 1/(N*dt) = 0.00095 * 100 Hz for simulation
	df = 1/(N_mu*dt) = 0.0019 * 100 Hz for mutest

*/
