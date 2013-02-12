//file "params.hh"

const double C_rate= .2; // firing rate [Hz]
const double C_dt= .01; // timestep [sec]
const int C_ndt=262144;//524288;//131072;//65536; // number of timesteps in simulation
const int C_ndt_mu=131072;//65536;//32768;//16384; // number of timesteps in bisection
const double C_eps=.2; // weight of Input-Current
const double C_tau_r=.1; // absolute refractory period [sec]
const int C_N_neuron=1000 ; // number of neurons per generation
const int C_N_Gen=5 ; // number of layers in the feedforward-networt (generations)

const int C_tau_r__dt=C_tau_r/C_dt; // absolute refractory period [#timesteps]
const double C_T_max=C_dt*(C_ndt-1); // simulation-time
const double C_T_mu=C_dt*(C_ndt_mu-1); // bisection simulation-time
const double C_tol_mu=1./C_T_mu; // tolerance for bisection
const int C_size_powspe=C_ndt/2+1; // size of vector for powerspectra

const double M2_PI=2*M_PI; // 2x pi
