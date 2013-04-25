// file "cl_sim.cc"

#define DEFINE_CLASS
#include "cl_sim.hh"
#undef DEFINE_CLASS

void ISI::lif_neuron(const double mu, const double* I_diff, const int ndt) // solve the differential equation tau_m * dv/dt = -v(t) + mu + eps*I(t) --> (tau_m=1) with fire-and-reset rule (v=1 --> spike at t --> v(t+tau_r)=0) and save the interspike-interval-times
{
	double v=0.;
	const double sum=mu+C_eps*C_rate, fak=C_eps*sqrt(C_rate)/C_N_neuron;
	int T=0;
	isi_.clear();
	for (unsigned int t=0; t<ndt; t++){	
		v+= (-v+sum+fak*I_diff[t])*C_dt; // Euler-step
		if (v>=1.){ // fire-and-reset rule
	      		isi_.push_back(T); // append 'T' to vector 'isi'
			T=C_tau_r__dt; // initial value for waiting time: absolute refractory period
			v=0.; // resetting the voltage
			t+=C_tau_r__dt; // nothing happens in absolute refractory period
	    	} else T++;
	}  
}

void ISI::calc_rho_k(double* rho_k, const unsigned int i_trial, const ISI& isi_train) // calculate the serial correlation coefficient rho_k until lag k_max-1
{
	double tmp1, tmp2, X1, X2, C, C1, C2;
	const unsigned int N=isi_train.isi_.size();

	for (unsigned int k = 0; k < C_rho_k_max; k++){	
		tmp1=N-k; tmp2=tmp1-5; X1=0; X2=0; C=0; C1=0; C2=0;
		#pragma omp parallel for \
			reduction(+: X1, X2, C, C1, C2)
		for (unsigned int i = 5; i < tmp1; i++){ // forget the first 5 ISIs
			X1+=isi_train.isi(i)*C_dt/tmp2;
			X2+=isi_train.isi(i+k)*C_dt/tmp2;
			C+=isi_train.isi(i)*C_dt*isi_train.isi(i+k)*C_dt/tmp2;
			C1+=(isi_train.isi(i)*C_dt)*(isi_train.isi(i)*C_dt)/tmp2;
			C2+=(isi_train.isi(i+k)*C_dt)*(isi_train.isi(i+k)*C_dt)/tmp2;
		}
		C-=X1*X2;
		C1-=X1*X1;
		C2-=X2*X2;
		rho_k[k+i_trial*C_rho_k_max]=C/sqrt(C1*C2);
	}		
}


