//functs.hh

#include "cl_sim.hh"
#include <fftw3.h>
#include <time.h>

void powerspectrum(double* spect, const ISI& isi_train); // calculate the powerspectrum of realisation "isi_train"

double mutest(const double* I_diff); // calculate mu for given rate r_0 (bisection)


void calc_I_diff(double* I_diff, double const* S, const fftw_plan p, const unsigned long int init); // calculate diffusion-current for one realisation

void dzeros(double *arr, const int n); // write zeros in double-array

void cpNcl_d_arr(double *arr_src, double *arr_trgt, const int size); // copy source-array to target-array and clear (set 0 values) source-array

double average(const double *arr, const int n); // calculate average of array elements

double variance(const double *arr, const int n, const double avg); // calculate variance of array elements

void whitenoise_PS(double *arr, const double r_0); // generate bandlimited white-noise powerspectrum with rate r_0 und cut-off frequency f_c=1/(2*C_dt)

void ornstein_PS(double *arr); // generate powerspectrum of an ohrnstein-uhlenbeck-process

double int_powspe(const double* powspe); // calculate the integral of the powerspectrum = variance^2

/***********************************************************************/

void safe_powspe(const double* powspe, const double* rho_avg, const double* rho_var, const double rate, const double sigma, const double CV, const double mu, const double rate_var, const double sigma_var, const double CV_var, const double mu_var, const int gen, const char* date); // saving powerspectrum to file for fix_r_cv=0 & fix_none=0

//void safe_powspe(const double* powspe, const double* rho_avg, const double* rho_var, const double rate, const double sigma, const double CV, const double rate_var, const double sigma_var, const double CV_var, const int gen, const char* date); // saving powerspectrum to file for fix_none=1

//void safe_powspe(const double* powspe, const double* rho_avg, const double* rho_var, const double rate, const double sigma, const double CV, const double mu, const double eps, const double rate_var, const double sigma_var, const double CV_var, const double mu_var, const double eps_var, const int gen, const char* date); // saving powerspectrum to file for fix_r_cv=1


/*double* mu_eps_test(double* mu_eps, const double r_0, const double cv_0, const int N_neuron, const double T_test, const double tol_mu, const double dt, const int tau_r__dt, const int N, const double* I_diff); // calculate mu and eps_diff for given rate r_0 and cv_0(bisection-algorithm)*/

