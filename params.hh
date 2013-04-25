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
