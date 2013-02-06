// file "main.cc"

#include "cl_sim.hh"
#include "cl_sim.cc"
//#include "powerspectrum.cc"
//#include "mutest.cc"

void dzeroes(double *arr, int n) // Nullen erzeugen f. double-Array
{
	for (unsigned int i=0;i<n;i++)
		arr[i]=0;
}

void whitenoise_PS(double *arr, int n, double r_0) // Erzeugt White-Noise Powerspectrum mit vorgegebener Rate
{
	for (unsigned int i=0;i<n;i++)
		arr[i]=r_0;
}

void cpNcl_d_arr(double *arr_src, double *arr_trgt, int size) // Copy Source-Array to Target-Array and clear (set 0 values) Source-Array
{
	for (unsigned int i = 0; i < size; i++)
	{
		arr_trgt[i]=arr_src[i];
		arr_src[i]=0;
	}
}

void calc_I_diff(double* I_diff, double const* S, int N,double dt, fftw_plan p, unsigned long int id, time_t now)
{
	int size=N/2;
	double df=1/(dt*N);
	double T_max= dt*(N-1);
	double a_n, phi_n;
	unsigned long int init= id*static_cast<unsigned long>(now);
	srand48(init);
  	gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set (rng,init);

	I_diff[0]=0;
	I_diff[size]=gsl_ran_gaussian_ziggurat(rng,sqrt(T_max*S[size]));
	for (unsigned int i=1;i<size;i++){
		a_n=gsl_ran_gaussian_ziggurat(rng,sqrt(T_max*S[i]));
		phi_n=M2_PI*drand48();
		I_diff[i]=a_n*cos(phi_n);
		I_diff[N-i]=a_n*sin(phi_n);
	}
	gsl_rng_free (rng);

	fftw_execute(p);
	
	double mean=0;
	double std_dev=0;

	for (unsigned int i = 0; i < N; i++)
	{
		I_diff[i]=df*I_diff[i];
		mean+=I_diff[i];
		std_dev+=I_diff[i]*I_diff[i];
	}

	cout << "Mean I= " << (mean/N) << "\t Std.dev. I= " << (std_dev/N) << endl;
}

double int_powspe(double* powspe, int size, double df)
{
	double integral=0;
	for (unsigned int i = 1; i < (size-1); i++)
	{
		integral+=powspe[i];
	}
	return (2*integral*df);
}

int main()
{
	double powspe_old[C_size_powspe]={0.0};
//		powspe_old=new double[C_size_powspe];
	double powspe_new[C_size_powspe]={0.0};
//		powspe_new=new double[C_size_powspe];
	double* S_temp;
	double* I_input;
	double mu;
	double mu_gen=0.;
	double rate_gen=0.;
	ISI interval(C_T_max, C_rate);
	
	unsigned long int rand_id=0;

	stringstream buffer;
	ofstream file;
	char date[18];
	const time_t now=time(NULL);
	strftime(date,18, "%Y-%m-%d_%H-%M",localtime(&now));
	stringstream filename;
	string filename_tmp;

/* create plan for fft in "calc_I_diff" */
	I_input=new double[C_ndt];
	fftw_plan plan_I_input = fftw_plan_r2r_1d(C_ndt, I_input, I_input, FFTW_HC2R, FFTW_MEASURE | FFTW_DESTROY_INPUT);

/* T */	double time, tstart;

	/* generate Start-Powerspectrum */
		whitenoise_PS(powspe_old,C_size_powspe,C_rate);
		dzeroes(powspe_new, C_size_powspe);

		for (unsigned int i_gen = 0; i_gen < C_N_Gen; i_gen++)
		{
/* T */			cout << "Gen: " << i_gen << endl;			
/* T */			tstart=clock();
			cout << "sigma^2 (integral S)= " << int_powspe(powspe_old, C_size_powspe, (1./(C_ndt*C_dt))) << endl;
			for (unsigned int i_neuron = 0; i_neuron < C_N_neuron; i_neuron++)
			{
/* T */				cout << "Neuron: " << i_neuron << endl;			
			/* generate Random-Numbers&Diffusion-Current*/
//				I_input=new double[C_ndt];
				rand_id++;
				calc_I_diff(I_input, powspe_old, C_ndt,C_dt, plan_I_input, rand_id, now);

//				interval.clearISI();
			/* define mu */
				mu=mutest(C_rate, C_T_mu, C_tol_mu, C_dt, C_taum, C_eps, C_ndt_mu, I_input);

			/* create ISI-train */
				interval.lif_neuron(mu, C_dt, C_taum, C_eps, C_ndt, I_input);
//				delete[] I_input;

			/* calculate Powerspectrum - normalized to the rate */
				S_temp= new double[C_size_powspe];
				powerspectrum(S_temp, interval, C_ndt, C_dt);

				for (unsigned int i_S = 0; i_S < C_size_powspe; i_S++)
				{
					powspe_new[i_S]+=S_temp[i_S]/C_N_neuron;
				}			
				delete[] S_temp;
			/* calculate averages of rate and mu */
				rate_gen+=interval.rate()/C_N_neuron;
				mu_gen+=mu/C_N_neuron;
			}

/* T */			time=clock()-tstart;
/* T */			cout << "Time im msec: " << (1000*time/CLOCKS_PER_SEC) << endl;


		/* safe Powerspectrum, mu's and rate's*/
//* T */			cout << "Start Saving " << endl;
//* T */			tstart = clock();
			buffer << "dt\t" << "N\t" << "tau\t" << "C_eps\t" << "C_N_neuron\n" << C_dt << "\t" << C_ndt << "\t" << C_taum << "\t" << C_eps << "\t" << C_N_neuron << "\n\nmu\t\trate\n" << mu_gen << "\t" << rate_gen << "\n\n";

/*			for (unsigned int i_safe = 0; i_safe < C_N_neuron; i_safe++)
			{
			buffer << "mu= " << mu_gen << "\t r_0= " << rate_gen << "\n \n";
			}
*/
			for (unsigned int i_safe = 0; i_safe < C_size_powspe; i_safe++)
			{
				buffer << (powspe_new[i_safe]/rate_gen) << endl;
			}

			filename << "data/" << date << "__" << i_gen << ".dat";
			filename_tmp=filename.str();

			file.open(filename_tmp.c_str());
				file << buffer.str();
			file.close();

			filename.str("");
			buffer.str("");

//* T */			time=clock()-tstart;
//* T */			cout << "Time in msec: " << (1000*time/CLOCKS_PER_SEC) << endl;

			cpNcl_d_arr(powspe_new, powspe_old, C_size_powspe);
			rate_gen=0.;
			mu_gen=0.;
		}
	
	delete[] I_input;
	fftw_destroy_plan(plan_I_input);
	fftw_forget_wisdom();

//	delete[] powspe_new;
//	delete[] powspe_old;

	return 0;
}
