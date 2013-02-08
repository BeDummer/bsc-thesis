//file "params"

const double C_rate= .5; // [Hz]
const double C_dt= .01; // [sec]
const int C_ndt=262144;//524288;//131072; // 65536;
const int C_ndt_mu=131072;//65536;//32768; // 16384 ;
const double C_D=0.11;
const double C_tau_r=.1;
const int C_N_neuron=100 ;
const int C_N_Gen=1 ;

const int C_tau_r__dt=C_tau_r/C_dt;
const double C_eps= sqrt(2.*C_D);
const double C_T_max=C_dt*(C_ndt-1);
const double C_T_mu=C_dt*(C_ndt_mu-1);
const double C_tol_mu=1./C_T_mu;
const int C_size_powspe=C_ndt/2+1;

const double M2_PI=2*M_PI;
