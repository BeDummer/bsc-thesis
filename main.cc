// file "main.cc"

#include "cl_sim.hh"
#include "cl_sim.cc"
/***********************************************************************/
void dzeros(double *arr, const int n) // write zeros in double-array
{
	for (unsigned int i=0;i<n;i++)
		arr[i]=0;
}
/***********************************************************************/
void cpNcl_d_arr(double *arr_src, double *arr_trgt, const int size) // copy source-array to target-array and clear (set 0 values) source-array
{
	for (unsigned int i = 0; i < size; i++)
	{
		arr_trgt[i]=arr_src[i];
		arr_src[i]=0;
	}
}
/***********************************************************************/
void whitenoise_PS(double *arr, const int n, const double r_0) // generate bandlimited white-noise powerspectrum with rate r_0 und cut-off frequency f_c=1/(2*C_dt)
{
	for (unsigned int i=0;i<n;i++)
		arr[i]=r_0;
}
/***********************************************************************/
void ornstein_PS(double *arr, const int N, const double df, const double tau, const double D) // generate powerspectrum of an ohrnstein-uhlenbeck-process
{
	double fak1=4.*tau*D, fak2=4.*M_PI*M_PI*df*df*tau*tau;
	for (unsigned int i = 0; i < N; i++)
	{
		arr[i]=fak1/(1+i*i*fak2);
	}
}
/***********************************************************************/
double int_powspe(const double* powspe, const int size, const double df) // calculate the integral of the powerspectrum = variance^2
{
	double integral=0;
	for (unsigned int i = 1; i < (size-1); i++)
	{
		integral+=powspe[i];
	}
	return (2*integral*df);
}
/***********************************************************************/
void safe_powspe(const double* powspe, const double rate, const double sigma, const double mu, const int gen, const char* date) // saving powerspectrum to file
{
	stringstream buffer, filename;
	string filename_tmp;
	ofstream file;

/* write Powerspectrum, mu, rate and constants into buffer */
	buffer << "dt\t" << "N\t" << "eps\t" << "tau_r\t" << "N_neuron\n" << C_dt << "\t" << C_ndt << "\t" << C_eps << "\t" << C_tau_r << "\t" << C_N_neuron << "\n\nmu\t\trate\t\tsigma\t\tCV\n" << mu << "\t" << rate << "\t" << sigma << "\t" << (sigma*rate) << "\n\n";

	for (unsigned int i = 0; i < C_size_powspe; i++)
	{
		buffer << (powspe[i]) << endl;
	}

	filename << "data/" << date << "__" << gen << ".dat";
	filename_tmp=filename.str();

/* saving data to file from buffer */
	file.open(filename_tmp.c_str());
		file << buffer.str();
	file.close();
}
/***********************************************************************/

int main()
{
/* variables for calculation */
	double powspe_old[C_size_powspe]={0.0};
	double powspe_new[C_size_powspe]={0.0};
	double* S_temp;
	double* I_input;
	double mu, mu_gen=0., rate_gen=0., sigma_gen=0.;
	ISI interval(C_T_max, C_rate);

/* variables for initialization of random number generators */	
	unsigned long int rand_id=1;
	const time_t now=time(NULL);

/* variables for saving data*/
	stringstream buffer;
	ofstream file;
	char date[18];
	strftime(date,18, "%Y-%m-%d_%H-%M",localtime(&now));
	stringstream filename;
	string filename_tmp;

/* create plan for fft in "calc_I_diff" */
	I_input=new double[C_ndt];
	fftw_plan plan_I_input = fftw_plan_r2r_1d(C_ndt, I_input, I_input, FFTW_HC2R, FFTW_MEASURE | FFTW_DESTROY_INPUT);

//* T */	double time, tstart;

	/* generate start-powerspectrum */
//		ornstein_PS(powspe_old, C_size_powspe, 1./(C_ndt*C_dt), C_tau, C_D);
		whitenoise_PS(powspe_old,C_size_powspe,C_rate);
		safe_powspe(powspe_old, C_rate, log(-1), log(-1), 0, date);
		dzeros(powspe_new, C_size_powspe);

		for (unsigned int i_gen = 1; i_gen <= C_N_Gen; i_gen++)
		{
/* T */			cout << "Gen: " << i_gen << endl;			
//* T */			tstart=clock();
//* T */			cout << "sigma^2 (integral S)= " << int_powspe(powspe_old, C_size_powspe, (1./(C_ndt*C_dt))) << endl;
			for (unsigned int i_neuron = 0; i_neuron < C_N_neuron; i_neuron++)
			{
/* T */				cout << "Neuron: " << i_neuron << endl;			
			/* generate Random-Numbers&Diffusion-Current*/
				rand_id++;
				calc_I_diff(I_input, powspe_old, C_ndt, C_dt, C_eps, C_rate, plan_I_input, rand_id, now);

			/* define mu */
				mu=mutest(C_rate, C_T_mu, C_tol_mu, C_dt, C_tau_r__dt, C_ndt_mu, I_input);

			/* create ISI-train */
				interval.lif_neuron(mu, C_dt, C_tau_r__dt, C_ndt, I_input);

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
				sigma_gen+=interval.std_dev(C_dt)/C_N_neuron;
				mu_gen+=mu/C_N_neuron;
			}

//* T */			time=clock()-tstart;
//* T */			cout << "Time im msec: " << (1000*time/CLOCKS_PER_SEC) << endl;

		/* saving powerspectrum, rate, sigma, mu and constants to file */
			safe_powspe(powspe_new, rate_gen, sigma_gen, mu_gen, i_gen, date);
		
		/* resetting variables */
			cpNcl_d_arr(powspe_new, powspe_old, C_size_powspe);
			rate_gen=0.;
			sigma_gen=0.;
			mu_gen=0.;
		}
	delete[] I_input;
	fftw_destroy_plan(plan_I_input);
	fftw_forget_wisdom();

	return 0;
}
