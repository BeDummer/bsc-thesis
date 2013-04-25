// functs.cc

#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "params.hh"
#include "cl_sim.hh"
	using namespace std;
#include <omp.h>

void powerspectrum(double* spect, const ISI& isi_train) // calculate the powerspectrum of realisation "isi_train"
{
	const unsigned int isi_size=isi_train.isi_.size();
	double train[C_ndt];
	const double fak=C_dt/(C_ndt-1), rec_dt=1/C_dt; 

	unsigned int count=0, j=0;
//	fftw_plan_with_nthreads(omp_get_max_threads());
	fftw_plan plan_fft;
	static char* fft_wisdom;
	static bool wisdom=false;

/* check for wisdom from earlier ffts*/
	if (wisdom){
		fftw_import_wisdom_from_string(fft_wisdom);
		plan_fft = fftw_plan_r2r_1d(C_ndt, train, train, FFTW_R2HC, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
	} else 	plan_fft = fftw_plan_r2r_1d(C_ndt, train, train, FFTW_R2HC, FFTW_MEASURE | FFTW_DESTROY_INPUT);

/* compute spiketrain from isi-times-train */
	for (unsigned int i=0; i<C_ndt; i++) {
		if (count<isi_size){
			if (j<isi_train.isi(count)) {
				train[i]=0;
				j++;
			} else {
				train[i]=rec_dt;
				count++;
				j=0;
			}
		} else 	train[i]=0;
	}

/* fourier-transform the spiketrain */
	fftw_execute(plan_fft);

	spect[0]=train[0]*train[0]*fak;
	spect[C_size_powspe-1]=train[C_ndt/2]*train[C_ndt/2]*fak;

	#pragma omp parallel for
	for (unsigned int i=1; i<(C_size_powspe-1); i++)
		spect[i]=(train[i]*train[i]+train[C_ndt-i]*train[C_ndt-i])*fak; // calculate the absolute value squared divided by the simulation time

/* save wisdom of the first fft */
	if (!wisdom) {
		fft_wisdom=fftw_export_wisdom_to_string();
		wisdom=true;
	}

	fftw_destroy_plan(plan_fft);
}

/**************************************************************/

double mutest(const double* I_diff) // calculate mu for given rate r_0 (bisection-algorithm)
{
	ISI test(C_T_mu);
	double min=-10. , max=5. , mid=(max+min)/2.;
	
	test.lif_neuron(mid,I_diff,C_ndt_mu);
	double r=test.rate();
	while (fabs(C_rate-r)>C_tol_mu) {
		(r>C_rate ? max : min) = mid;
		mid=(max+min)/2;
		test.lif_neuron(mid,I_diff,C_ndt_mu);
		r=test.rate();
	}
/* T */	cout << mid << "\t" << r << endl;
	return mid;
}

/**************************************************************/

void calc_I_diff(double* I_diff, double const* S, const fftw_plan p, const unsigned long int init) // calculate diffusion-current for one realisation
{
	const int size=C_ndt/2;
	const double fak2= C_T_max*C_eps*C_eps/C_N_neuron, fak1=4.*M_PI*M_PI*C_df*C_df*C_tau_s*C_tau_s; 
	double a_n, phi_n;
	unsigned int i;
	srand48(init);
  	gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set (rng,init);

/* generate current in fourier-domain with the statistics of the input-powerspectrum */
	I_diff[0]=0;
	I_diff[size]=gsl_ran_gaussian_ziggurat(rng,sqrt(fak2*S[size]/(1.+fak1*size*size)));
	#pragma omp parallel for \
		shared(I_diff, rng, fak2, S, fak1, M2_PI, C_ndt, size) private(i, a_n, phi_n)
	for (i=1;i<size;i++){
		#pragma omp critical
		{
			a_n=gsl_ran_gaussian_ziggurat(rng,sqrt(fak2*S[i]/(1.+fak1*i*i)));
			phi_n=M2_PI*drand48();
		}
		I_diff[i]=a_n*cos(phi_n);
		I_diff[C_ndt-i]=a_n*sin(phi_n);
	}
	gsl_rng_free(rng);

/* fourier-transform the current */
	fftw_execute(p);

/* multiply with df */
	#pragma omp parallel for
	for (i = 0; i < C_ndt; i++) I_diff[i]=C_df*I_diff[i];
}
/***********************************************************************/
void dzeros(double *arr, const int n) // write zeros in double-array
{
	#pragma omp parallel for
	for (unsigned int i=0;i<n;i++) arr[i]=0;
}
/***********************************************************************/
void cpNcl_d_arr(double *arr_src, double *arr_trgt, const int size) // copy source-array to target-array and clear (set 0 values) source-array
{
	#pragma omp parallel for
	for (unsigned int i = 0; i < size; i++) {
		arr_trgt[i]=arr_src[i];
		arr_src[i]=0;
	}
}
/***********************************************************************/
double average(const double *arr, const int n) // calculate average of array elements
{
	double tmp=0;
	#pragma omp parallel for \
		reduction(+: tmp)
	for (unsigned int i=0;i<n;i++) tmp +=arr[i]/n;
	return tmp;
}
/***********************************************************************/
double variance(const double *arr, const int n, const double avg) // calculate variance of array elements
{
	double tmp=0;
	#pragma omp parallel for \
		reduction(+: tmp)
	for (unsigned int i=0;i<n;i++) tmp +=(avg-arr[i])*(avg-arr[i])/(n-1);
	return tmp;
}
/***********************************************************************/
void whitenoise_PS(double *arr, const double r_0) // generate bandlimited white-noise powerspectrum with rate r_0 und cut-off frequency f_c=1/(2*C_dt)
{
	#pragma omp parallel for
	for (unsigned int i=0;i<C_size_powspe;i++) arr[i]=r_0;
}
/***********************************************************************/
void ornstein_PS(double *arr) // generate powerspectrum of an ohrnstein-uhlenbeck-process
{
	double fak1=4.*C_tau_ou*C_D_ou, fak2=4.*M_PI*M_PI*C_df*C_df*C_tau_ou*C_tau_ou;
	#pragma omp parallel for
	for (unsigned int i = 0; i < C_size_powspe; i++) arr[i]=fak1/(1+i*i*fak2);
}
/***********************************************************************/
double int_powspe(const double* powspe) // calculate the integral of the powerspectrum = variance^2
{
	double integral=0;
	#pragma omp parallel for \
		reduction(+: integral)
	for (unsigned int i = 1; i < (C_size_powspe-1); i++) integral+=powspe[i];

	return (2*integral*C_df);
}
/***********************************************************************/
void safe_powspe(const double* powspe, const double* rho_avg, const double* rho_var, const double rate, const double sigma, const double CV, const double mu, const double rate_var, const double sigma_var, const double CV_var, const double mu_var, const int gen, const char* date) // saving powerspectrum to file for fix_r_cv=0 & fix_none=0
{
	stringstream buffer, filename;
	string filename_tmp;
	ofstream file;
	unsigned int i;

// write Powerspectrum, mu, rate and constants into buffer
	buffer << "dt\t" << "N\t" << "eps\t" << "tau_r\t" << "tau_s\t" << "N_neuron\t" << "N_trials\t" << "k_max\t" << "r_0\n"
		<< C_dt << "\t" << C_ndt << "\t" << C_eps << "\t" << C_tau_r << "\t" << C_tau_s << "\t" << C_N_neuron << "\t" << C_N_trials << "\t" << C_rho_k_max << "\t" << C_rate
		<< "\n\n mu\t\trate\t\tsigma\t\tCV\n"
		<< mu << "\t" << rate << "\t" << sigma << "\t" << CV << "\n" << mu_var << "\t" << rate_var << "\t" << sigma_var << "\t" << CV_var << "\n\n rho_avg \t rho_var\n";

	for (i = 0; i < C_rho_k_max; i++) buffer << rho_avg[i] << "\t" << rho_var[i] << endl;

	buffer << endl;

	for (i = 0; i < C_size_powspe; i++) buffer << (powspe[i]) << endl;

	filename << "data/" << date << "_id-1_r-" << C_rate << "_e-" << C_eps << "__" << gen << ".dat";
	filename_tmp=filename.str();

// saving data to file from buffer 
	file.open(filename_tmp.c_str());
		file << buffer.str();
	file.close();
}
