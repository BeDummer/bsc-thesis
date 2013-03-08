// file "cl_sim.cc"

#include "cl_sim.hh"

/*************************************************/

void ISI::lif_neuron(const double mu, const double r_0, const double eps, const int N_neuron, const double dt, const int tau_r__dt, const int N_step, const double* I_diff) // solve the differential equation tau_m * dv/dt = -v(t) + mu + eps*I(t) --> (tau_m=1) with fire-and-reset rule (v=1 --> spike at t --> v(t+tau_r)=0) and save the interspike-interval-times
{
	double v=0.;
	const double sum=mu+eps*r_0, fak=eps*sqrt(r_0)/N_neuron;
	int T=0;
	isi_.clear();
	for (unsigned int t=0; t<N_step; t++)
		{
			v+= (-v+sum+fak*I_diff[t]))*dt; // Euler-step
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

void ISI::calc_rho_k(const unsigned int k_max, const double dt, double* rho_k, const unsigned int neuron, const ISI& isi_train)
{
	double tmp, X1, X2, C, C1, C2;
	const unsigned int N=isi_train.isi_.size();

for (unsigned int k = 0; k < k_max; k++)
{
	tmp=N-k; X1=0; X2=0; C=0; C1=0; C2=0;
	for (unsigned int i = 0; i < tmp; i += 1)
	{
		X1+=isi_train.isi(i)*dt/tmp;
		X2+=isi_train.isi(i+k)*dt/tmp;
		C+=isi_train.isi(i)*dt*isi_train.isi(i+k)*dt/tmp;
		C1+=pow((isi_train.isi(i)*dt),2)/tmp;
		C2+=pow((isi_train.isi(i+k)*dt),2)/tmp;
	}
	C-=X1*X2;
	C1-=pow(X1,2);
	C2-=pow(X2,2);
	rho_k[k+neuron*k_max]=C/sqrt(C1*C2);
}

/*	for (unsigned int k = 0; k < k_max; k++) NON-ESTIMATOR FOR RHO_K
	{
		tmp=0;
		for (unsigned int i = 0; (i+k) < isi_size; i++)
		{
			tmp+=(isi_train.isi(i)*dt)*(isi_train.isi(i+k)*dt)/(isi_size-k);
		}
		rho_k[k+neuron*k_max]=(tmp-T_sqr)/variance;
	}
*/	
}

/*****************************************************

inline double window(int i, int size) // Welch window function
{
	double x=1-pow((2*i/size-1),2);
	return x;
}

/*****************************************************/

void powerspectrum(double* spect, const ISI& isi_train, const int N, const double dt) // calculate the powerspectrum of realisation "isi_train"
{
	const unsigned int size_powspe=N/2+1, isi_size=isi_train.isi_.size();
	double* train = new double[N];
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
	delete[] train;
}

/**************************************************************/

double mutest(const double r_0, const double eps, const int N_neuron const int tau_r__dt, const double T_test, const double tol_mu, const double dt, const int N_step, const double* I_diff) // calculate mu for given rate r_0 (bisection-algorithm)
{
	ISI test(T_test, r_0);
	double min=-2. , max=2. , mid=(max+min)/2.;
	
	test.lif_neuron(mid,r_0,eps,dt,tau_r__dt,N,I_diff);
	double r=test.rate();
	while (fabs(r_0-r)>tol_mu) {
		(r>r_0 ? max : min) = mid;
		mid=(max+min)/2;
		test.lif_neuron(mid,r_0,eps,N_neuron,dt,tau_r__dt,N_step,I_diff);
		r=test.rate();
//* T */		cout << mid << "\t" << r << endl;
	}
/* T */	cout << mid << "\t" << r << endl;
	return mid;
}

/**************************************************************/

double* mu_eps_test(const double r_0, const double cv_0, const int N_neuron, const int tau_r__dt, const double T_test, const double tol_mu, const double dt, const int N_step, const double* I_diff) // calculate mu and eps_diff for given rate r_0 and cv_0(bisection-algorithm)
{
	ISI test(T_test, r_0);
	const double tol=pow(10.,-6.);
	double tol_cv=5*tol_mu;
	int count=0;
	double min[2], max[2], mid[2], *mu_eps;
	min[1]=-2.;
	min[2]=0.;
	max[1]=2.;
	max[2]=1.;
	mid[1]=(max[1]+min[1])/2.;
	mid[2]=(max[2]+min[2])/2.;

	test.lif_neuron(mid[1],r_0,mid[2],N_neuron,dt,tau_r__dt,N_step,I_diff);
	double r=test.rate();
	double cv=r*sqrt(test.var(dt));
	while (fabs(cv_0-cv)>tol_cv || fabs(r_0-r)>tol_mu) {
		count++;
		if (fabs(max[2]-min[2])<tol)
		{
			min[2]/=10.;
			max[2]*=10.;
			mid[2]=(max[2]+min[2])/2.;
			test.lif_neuron(mid[1],r_0,mid[2],N_neuron,dt,tau_r__dt,N_step,I_diff);
			r=test.rate();
		}
		while (fabs(r_0-r)>tol_mu) {
			(r>r_0 ? max[1] : min[1]) = mid[1];
			mid[1]=(max[1]+min[1])/2.;
			test.lif_neuron(mid[1],r_0,mid[2],N_neuron,dt,tau_r__dt,N_step,I_diff);
			r=test.rate();
		}
		cv=r*sqrt(test.var(dt));
		(cv>cv_0 ? max[2] : min[2]) = mid[2];
		mid[2]=(max[2]+min[2])/2.;
		test.lif_neuron(mid[1],r_0,mid[2],N_neuron,dt,tau_r__dt,N_step,I_diff);
		r=test.rate();		
		cv=r*sqrt(test.var(dt));
		min[1]-=.5;
		max[1]+=.5;

		if (count>100) {
			tol_cv*=2.;
			if (count>1000) {
/* T */				cout << "!!!No Convergence in mu_eps_test!!!" << endl;
				mid[1]=log(-1);
				mid[2]=log(-1);
				mu_eps=mid;
				break;				
			}
		}	
//* T */		cout << mid[1] << "\t" << r << "\t" << mid[2] << "\t" << cv << endl;
	}
/* T */	cout << mid[1] << "\t" << r << "\t" << mid[2] << "\t" << cv << endl;
		
	mu_eps=mid;
	return mu_eps;
}

/**************************************************************/

void calc_I_diff(double* I_diff, double const* S, const double eps, const int N_neuron, const int N_step, const double dt, const double r_0, const double tau_s, const fftw_plan p, const unsigned long int id, const time_t now) // calculate diffusion-current for one realisation
{
	const int size=N_step/2;
	const double df=1/(dt*N_step), fak2= dt*(N_step-1)*pow(eps,2)/N_neuron, fak1=4*MSQR_PI*pow(df,2)*pow(tau_s,2); 
	double a_n, phi_n;
	unsigned long int init= id+static_cast<unsigned long>(now);
	srand48(init);
  	gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set (rng,init);

/* generate current in fourier-domain with the statistics of the input-powerspectrum */
	I_diff[0]=0;
	I_diff[size]=gsl_ran_gaussian_ziggurat(rng,sqrt(fak2*S[size]/(1+fak1*pow(size,2))));
	for (unsigned int i=1;i<size;i++){
		a_n=gsl_ran_gaussian_ziggurat(rng,sqrt(fak2*S[i]/(1+fak1*pow(i,2))));
		phi_n=M2_PI*drand48();
		I_diff[i]=a_n*cos(phi_n);
		I_diff[N_step-i]=a_n*sin(phi_n);
	}
	gsl_rng_free (rng);

/* fourier-transform the current */
	fftw_execute(p);
	
//* T */	double mean=0;
//* T */	double std_dev=0;

/* multiply with df */
	for (unsigned int i = 0; i < N_step; i++)
	{
		I_diff[i]=df*I_diff[i];
//* T */		mean+=I_diff[i];
//* T */		std_dev+=I_diff[i]*I_diff[i];
	}
//* T */	cout << "Mean I= " << (mean/N_step) << "\t Std.dev. I= " << (std_dev/N_step) << endl;
}

