// file "cl_sim.hh"

#ifndef _cl_sim_hh_
#define _cl_sim_hh_

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "params.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

/*******************************************************************/

class ISI {
	private:
		vector<int> isi_;
		double T_max_;
	public:
		ISI(double T_max, double r_0):
			T_max_(T_max)
			{isi_.reserve(static_cast<int>(T_max*r_0*1.5));};
		int isi(int i) const {return isi_[i];};
		void lif_neuron(const double mu, const double dt, const double eps, const int tau_r__dt, const int N, const double* I_diff);
		double rate();
	friend void powerspectrum(double* , const ISI&, const int, const double);
};

void powerspectrum(double* spect, const ISI& isi_train, const int N, const double dt); // calculate the powerspectrum of realisation "isi_train"

double mutest(const double r_0, const double T_test, const double tol_mu, const double dt, const double eps, const int tau_r__dt, const int N, const double* I_diff); // calculate mu for given rate r_0 (bisection)

void calc_I_diff(double* I_diff, double const* S, const int N, const double dt, const fftw_plan p, const unsigned long int id, const time_t now); // calculate diffusion-current for one realisation
/**************************************************************************/

#endif //_cl_sim_hh_
