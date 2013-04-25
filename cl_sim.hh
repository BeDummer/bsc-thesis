// file "cl_sim.hh"

#ifdef DEFINE_CLASS
	#include "params.hh"
#endif
#include <math.h>
#include <vector>
	using namespace std;

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
	for (unsigned int i = 0; i < size; i++){
		tmp+=(T-isi_[i]*C_dt)*(T-isi_[i]*C_dt)/(size-1);
	}
	return tmp;
}
