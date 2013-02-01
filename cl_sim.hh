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
#include "params"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

/*
class Spike {
	private:
		short unsigned int* train_;
		int ndt_;
	public:
		Spike(int N):
			ndt_(N)
			{train_ = new short unsigned int[N];};
		~Spike() {delete train_[];};
		void isi_spike(const ISI& isi_train);
		short unsigned int train(int i) const {return train_[i];};
};

class Power {
	private:
		double* spect_;
		short unsigned int* train_;
		int ndt_;
		void isi_spike(const ISI& isi_train);
	public:
		Power(int N, double init=0):
			ndt_(N)
			{spect_= new double[(N/2+1)]; for (unsigned int i=0;i<(N/2+1);i++) {spect_[i]=init;} }
		~Power() {delete spect_[];}
		double spect(int i) const {return spect_[i];};
		Power powerspectrum(int ndt, const Spike& train);
};
*/
/*************************************************************
class Random {
	private:
		double* a_n_; 
		double* phi_n_;
	public:
		Random(int N, double const* S);
		double a_n(int i) const {return a_n_[i];};
		double phi_n(int i) const {return phi_n_[i];};
		~Random();
};
*******************************************************************/

class ISI {
	private:
		vector<int> isi_;
		double T_max_;
	public:
		ISI(double T_max, double r_0):
			T_max_(T_max)
			{isi_.reserve(static_cast<int>(T_max*r_0*1.5));};
//		void clearISI() {isi_.clear();};
		int isi(int i) const {return isi_[i];};
		void lif_neuron(double mu, double dt, double taum, double eps, double N, const double* I_diff);
		double rate();
	friend void powerspectrum(double* , const ISI&, int, double);
};

void powerspectrum(double* spect, const ISI& isi_train, int N, double dt);

double mutest(double r_0, double T_test, double tol_mu, double dt, double taum, double eps, double N, const double* I_diff);

/**************************************************************************
inline Random::Random(int N, double const* S)
{
	int size=N/2;
	a_n_=new double[size];
	phi_n_=new double[size];
	srand48((long)time(NULL));
  	gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set (rng,(long)time(NULL));

	for (unsigned int i=0;i<size;i++){
		a_n_[i]=gsl_ran_gaussian_ziggurat(rng,sqrt(S[i+1]));
		phi_n_[i]=2*M_PI*drand48();
	}
	gsl_rng_free (rng);

}

inline Random::~Random()
{
	delete[] a_n_;
	delete[] phi_n_;
}
***************************************************************************/

#endif //_cl_sim_hh_
