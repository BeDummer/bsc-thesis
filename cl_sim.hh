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
		int isi(int i)const {return isi_[i];};
		double rate()const {return (isi_.size()/T_max_);};
		double CV(const double dt)const {return (rate()*var(dt));};
		double var(const double dt)const;
		void calc_rho_k(const unsigned int, const double, double*, const int, const ISI&);
		void lif_neuron(const double, const double, const double, const double , const int, const int, const double*);
	friend void powerspectrum(double* , const ISI&, const int, const double);
};

inline double ISI::var(const double dt)const
{
	double T_sqr=1./pow(ISI::rate(),2), tmp=0;
	int size=isi_.size();
	for (unsigned int i = 0; i < size; i++)
	{
		tmp+=pow((isi_[i]*dt),2)/size;
	}
	return (tmp-T_sqr);
}

void powerspectrum(double* spect, const ISI& isi_train, const int N, const double dt); // calculate the powerspectrum of realisation "isi_train"

double mutest(const double r_0, const double eps_avg, const double eps_diff, const int tau_r__dt, const double T_test, const double tol_mu, const double dt, const int N, const double* I_diff); // calculate mu for given rate r_0 (bisection)

double* mu_eps_test(const double r_0, const double cv_0, const double eps_avg, const double T_test, const double tol_mu, const double dt, const int tau_r__dt, const int N, const double* I_diff); // calculate mu and eps_diff for given rate r_0 and cv_0(bisection-algorithm)

void calc_I_diff(double* I_diff, double const* S, const int N, const double dt, const double r_0, const double tau_s, const fftw_plan p, const unsigned long int id, const time_t now); // calculate diffusion-current for one realisation
/**************************************************************************/

#endif //_cl_sim_hh_
