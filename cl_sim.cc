// file "cl_sim.cc"

#include "cl_sim.hh"

/*************************************************/

void ISI::lif_neuron(const double mu, const double dt, const int tau_r__dt, const int N, const double* I_diff) // solve the differential equation tau_m * dv/dt = -v(t) + mu + eps*I(t) --> (tau_m=1) with fire-and-reset rule (v=1 --> spike at t --> v(t+tau_r)=0) and save the interspike-interval-times
{
	double v=0.;
	int T=0;
	isi_.clear();
	for (unsigned int t=0; t<N; t++)
		{
			v+= (-v+mu+I_diff[t])*dt; // Euler-step
			if (v>=1.) // fire-and-reset rule
		    	{
		      		isi_.push_back(T); // append 'T' to vector 'isi'
				T=tau_r__dt; // initial value for waiting time: absolute refractory period
				v=0.; // resetting the voltage
				t+=tau_r__dt; // nothing happens in absolute refractory period
		    	} else
				T++;
		}  
}

double ISI::rate()
{
	double r=isi_.size()/T_max_;
	return r;
}

/******************************************************/

void powerspectrum(double* spect, const ISI& isi_train, const int N, const double dt) // calculate the powerspectrum of realisation "isi_train"
{
	const unsigned int size_powspe=N/2+1, isi_size=isi_train.isi_.size();
	double train[N];
	const double fak=dt/(N-1), rec_dt=1/dt;

	unsigned int count=0, j=0;

	fftw_plan plan_fft;
	static char* fft_wisdom;
	static bool wisdom=false;

/* check for wisdom from earlier ffts*/
	if (wisdom)
	{
		fftw_import_wisdom_from_string(fft_wisdom);
		plan_fft = fftw_plan_r2r_1d(N, train, train, FFTW_R2HC, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
	} else
		plan_fft = fftw_plan_r2r_1d(N, train, train, FFTW_R2HC, FFTW_MEASURE | FFTW_DESTROY_INPUT);

/* compute spiketrain from isi-times-train */
	for (unsigned int i=0; i<N; i++) {
		if (count<isi_size){
			if (j<isi_train.isi(count)) {
				train[i]=0;
				j++;
			} else {
				train[i]=rec_dt;
				count++;
				j=0;
			}
		} else 
			train[i]=0;
	}
//* T */		double rate=0;
//**/		for (unsigned int i = 0; i < N; i++)
//**/		{
//**/			rate+=train[i]/(N-1);
//**/		}
//* T */		cout << "Train-rate= " << rate << endl;

/* fourier-transform the spiketrain */
	fftw_execute(plan_fft);
	spect[0]=train[0]*train[0]*fak;
	spect[size_powspe-1]=train[N/2]*train[N/2]*fak;

	for (unsigned int i=1; i<(size_powspe-1); i++)
		spect[i]=(train[i]*train[i]+train[N-i]*train[N-i])*fak; // calculate the absolute value squared divided by the simulation time

/* save wisdom of the first fft */
	if (!wisdom) 
	{
		fft_wisdom=fftw_export_wisdom_to_string();
		wisdom=true;
	}

	fftw_destroy_plan(plan_fft);
}

/**************************************************************/

double mutest(const double r_0, const double T_test, const double tol_mu, const double dt, const int tau_r__dt, const int N, const double* I_diff) // calculate mu for given rate r_0 (bisection-algorithm)
{
	ISI test(T_test, r_0);
	double min=-5. , max=5. , mid=(max+min)/2.;
	
	test.lif_neuron(mid,dt,tau_r__dt,N,I_diff);
	double r=test.rate();
	while (fabs(r_0-r)>tol_mu) {
		if (r>r_0)
      			max=mid;
    		else
      			min=mid;
		mid=(max+min)/2;
		test.lif_neuron(mid,dt,tau_r__dt,N,I_diff);
		r=test.rate();
//* T */		cout << mid << "\t" << r << endl;
	}
/* T */	cout << mid << "\t" << r << endl;
	return mid;
}

/**************************************************************/

void calc_I_diff(double* I_diff, double const* S, const int N, const double dt, const double eps, const double r_0, const fftw_plan p, const unsigned long int id, const time_t now) // calculate diffusion-current for one realisation
{
	const int size=N/2;
	const double df=1/(dt*N), T_max= dt*(N-1); 
	double a_n, phi_n;
	unsigned long int init= id+static_cast<unsigned long>(now);
	srand48(init);
  	gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set (rng,init);

/* generate current in fourier-domain with the statistics of the input-powerspectrum */
	I_diff[0]=0;
	I_diff[size]=gsl_ran_gaussian_ziggurat(rng,sqrt(T_max*eps*eps*S[size]));
	for (unsigned int i=1;i<size;i++){
		a_n=gsl_ran_gaussian_ziggurat(rng,sqrt(T_max*eps*eps*S[i]));
		phi_n=M2_PI*drand48();
		I_diff[i]=a_n*cos(phi_n);
		I_diff[N-i]=a_n*sin(phi_n);
	}
	gsl_rng_free (rng);

/* fourier-transform the current */
	fftw_execute(p);
	
//* T */	double mean=0;
//* T */	double std_dev=0;

/* adding the average current */
	for (unsigned int i = 0; i < N; i++)
	{
		I_diff[i]=eps*r_0+df*I_diff[i];
//* T */		mean+=I_diff[i];
//* T */		std_dev+=I_diff[i]*I_diff[i];
	}
//* T */	cout << "Mean I= " << (mean/N) << "\t Std.dev. I= " << (std_dev/N) << endl;
}
