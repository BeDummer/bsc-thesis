// file "cl_sim.cc"

#include "cl_sim.hh"

/*************************************************/

void ISI::lif_neuron(double mu, double dt, double taum, double eps, double N, const double* I_diff)
{
	double cons=dt/taum, v=0.;
	int T=0;
	isi_.clear();
	for (unsigned int t=0; t<N; t++)
		{
			v+= (-v+mu+eps*I_diff[t])*cons; //Euler-step
			if (v>=1.)
		    	{
		      		isi_.push_back(T); // append 'T' to vector 'isi'
				T=0;
				v=0.;
		    	} else
				T++;
		}  
}

double ISI::rate()
{
	double r=isi_.size()/T_max_;
	return r;
}

/**********************************************************************
double ISI::I_diff(double time, double dt, int N, const Random& randval)
{
	double I=0;
	double df=1/(dt*N);
	for (unsigned int i=0;i<(N/2);i++){
		I+=randval.a_n(i) *cos(df*(i+1)*time-randval.phi_n(i));
	}

	return I;
}
*********************************************************************/


/******************************************************/


void powerspectrum(double* spect, const ISI& isi_train, int N, double dt)
{
	int size_powspe=N/2+1;
	double train[N];
	double fak=dt/(N-1), rec_dt=1/dt;

	unsigned int count=0, j=0, isi_size=isi_train.isi_.size();

	fftw_plan plan_fft;
	static char* fft_wisdom;
	static bool wisdom=false;

	if (wisdom)
	{
		fftw_import_wisdom_from_string(fft_wisdom);
		plan_fft = fftw_plan_r2r_1d(N, train, train, FFTW_R2HC, FFTW_WISDOM_ONLY | FFTW_DESTROY_INPUT);
	} else
		plan_fft = fftw_plan_r2r_1d(N, train, train, FFTW_R2HC, FFTW_MEASURE | FFTW_DESTROY_INPUT);

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

/* T */		double rate=0;
/**/		for (unsigned int i = 0; i < N; i++)
/**/		{
/**/			rate+=train[i]/(N-1);
/**/		}
/* T */		cout << "Train-rate= " << rate << endl;

	fftw_execute(plan_fft);
	spect[0]=train[0];
	spect[size_powspe-1]=train[N/2];
	for (unsigned int i=1; i<(size_powspe-1); i++)
		spect[i]=(train[i]*train[i]+train[N-i]*train[N-i])*fak;
	
	if (!wisdom) 
	{
		fft_wisdom=fftw_export_wisdom_to_string();
		wisdom=true;
	}

	fftw_destroy_plan(plan_fft);
}

/**************************************************************/

double mutest(double r_0, double T_test, double tol_mu, double dt, double taum, double eps, double N, const double* I_diff)
{
	ISI test(T_test, r_0);
	double min=-5. , max=5. , mid=(max+min)/2.;
	
	test.lif_neuron(mid,dt,taum,eps,N,I_diff);
	double r=test.rate();
	while (fabs(r_0-r)>tol_mu) {
		if (r>r_0)
      			max=mid;
    		else
      			min=mid;
		mid=(max+min)/2;
//		test.clearISI();
		test.lif_neuron(mid,dt,taum,eps,N,I_diff);
		r=test.rate();
/* T */		cout << mid << "\t" << r << endl;
	}
//* T */	cout << mid << "\t" << r << endl;
	return mid;
}

/*****************************************************

Power Power::powerspectrum(const Spike& spiketrain)
{
	Power S();
	int size_powspe=ndt_/2+1;
	fftw_plan p;
	fftw_complex *outfft;
	outfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_powspe);
	p = fftw_plan_dft_r2c_1d(ndt_, spiketrain.train, outfft, FFTW_ESTIMATE);

	fftw_execute(p);

	for (unsigned int i=0; i<size_powspe; i++)
		S.spect_[i]=outfft[i][0]*outfft[i][0]+outfft[i][1]*outfft[i][1];
	
	fftw_free(outfft);
	fftw_destroy_plan(p);
	
	return S;
	~S;
}

/*********************************************************

void Spike::isi_spike(const ISI& isi_train)
{
	unsigned int count=0, j=0, isi_size=isi_train.isi_.size();
	for (unsigned int i=0; i<ndt_; i++) {
		if (count<isi_size){
			if (j<isi_train.isi(count)) {
				train_[i]=0;
				j++;
			} else {
				train_[i]=1;
				count++;
				j=0;
			}
		} else 
			train_[i]=0;
}

/************************************************************/
