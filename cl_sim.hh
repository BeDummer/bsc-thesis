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
		void lif_neuron(const double mu, const double dt, const double taum, const double eps, const int N, const double* I_diff);
		double rate();
	friend void powerspectrum(double* , const ISI&, const int, const double);
};

void powerspectrum(double* spect, const ISI& isi_train, const int N, const double dt);

double mutest(const double r_0, const double T_test, const double tol_mu, const double dt, const double taum, const double eps, const int N, const double* I_diff);

/**************************************************************************/

#endif //_cl_sim_hh_
