// file "cl_sim.hh"

#ifdef DEFINE_CLASS
	#include "params.hh"
#endif
#include <math.h>
#include <vector>
	using namespace std;
/*******************************************************************/

class ISI {
	private:
		vector<int> isi_;
		double T_max_;
	public:
		ISI(double T_max):
			T_max_(T_max)
			{isi_.reserve(static_cast<int>(T_max*C_rate*1.5));};
		~ISI() {isi_.clear();};
		int isi(int i)const {return isi_[i];};
		double rate()const {return (isi_.size()/T_max_);};
		double CV()const {return (rate()*sqrt(var()));};
		double var()const;
		void calc_rho_k(double*, const unsigned int, const ISI&);
		void lif_neuron(const double, const double*, const int);
	friend void powerspectrum(double* , const ISI&);
};

inline double ISI::var()const
{
	double T=1./ISI::rate(), tmp=0.;
	int size=isi_.size();
	for (unsigned int i = 0; i < size; i++)
	{
		tmp+=(T-isi_[i]*C_dt)*(T-isi_[i]*C_dt)/(size-1);
	}
	return tmp;
}


/**************************************************************************
void ISI::lif_neuron(const double mu, const double r_0, const double eps, const int N_neuron, const double dt, const int tau_r__dt, const int N_step, const double* I_diff); // solve the differential equation tau_m * dv/dt = -v(t) + mu + eps*I(t) --> (tau_m=1) with fire-and-reset rule (v=1 --> spike at t --> v(t+tau_r)=0) and save the interspike-interval-times

void ISI::calc_rho_k(const unsigned int k_max, const double dt, double* rho_k, const unsigned int neuron, const ISI& isi_train); // calculate the serial correlation coefficient rho_k until lag k_max-1
*/
