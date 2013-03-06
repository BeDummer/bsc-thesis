// file "main.cc"

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
double average(const double *arr, const int n) // calculate average of array elements
{
	double tmp=0;
	for (unsigned int i=0;i<n;i++)
		tmp +=arr[i]/n;
	return tmp;
}
/***********************************************************************/
double variance(const double *arr, const int n, const double avg) // calculate variance of array elements
{
	double tmp=0;
	for (unsigned int i=0;i<n;i++)
		tmp +=pow(arr[i],2)/n;
	return (tmp-avg);
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
	double fak1=4.*tau*D, fak2=4.*MSQR_PI*pow(df,2)*pow(tau,2);
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
void safe_powspe(const double* powspe, const double* rho_avg, const double* rho_var, const double rate, const double sigma, const double mu, const double rate_var, const double sigma_var, const double mu_var, const int gen, const char* date) // saving powerspectrum to file
{
	stringstream buffer, filename;
	string filename_tmp;
	ofstream file;
	unsigned int i;

/* write Powerspectrum, mu, rate and constants into buffer */
	buffer << "dt\t" << "N\t" << "eps\t" << "tau_r\t" << "tau_s\t" << "N_neuron\t" << "k_max\t" << "r_0\t" << "CV_0\n"
		<< C_dt << "\t" << C_ndt << "\t" << C_eps << "\t" << C_tau_r << "\t" << C_tau_s << "\t" << C_N_neuron << "\t" << C_rho_k_max << "\t" << C_rate << "\t" << C_cv
		<< "\n\n mu\t\trate\t\tsigma\t\tCV\n"
		<< mu << "\t" << rate << "\t" << sigma << "\t" << (sqrt(sigma)*rate) << "\n" << mu_var << "\t" << rate_var << "\t" << sigma_var << "\t" << sqrt(pow((rate*sigma_var),2)/sigma+sigma*pow(rate_var,2)) << "\n\n rho_avg \t rho_var\n";

	for (i = 0; i < C_rho_k_max; i++)
	{
		buffer << rho_avg[i] << "\t" << rho_var[i] << endl;
	}

	buffer << "\n";

	for (i = 0; i < C_size_powspe; i++)
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

int main(int argc, char *argv[])
{
	bool arg_mu_cv = (argc>1 ? 1 : 0);

/* variables and constants for calculation */
	double *powspe_old, *powspe_new, *S_temp, *I_input;
	double *mu_gen, *rate_gen, *sigma_gen, *eps_gen, *rho_k_gen, *mu_eps, mu;
	double mu_avg, rate_avg, sigma_avg, eps_avg, rho_avg[C_rho_k_max]={0.0}, rho_var[C_rho_k_max]={0.0};

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
		powspe_old=new double[C_size_powspe];
//		ornstein_PS(powspe_old, C_size_powspe, 1./(C_ndt*C_dt), C_tau, C_D); safe_powspe(powspe_old, rho_avg, rho_var, C_rate, C_tau, C_D, log(-1), log(-1), log(-1), 0, date);
		whitenoise_PS(powspe_old,C_size_powspe,C_rate); safe_powspe(powspe_old, rho_avg, rho_var, C_rate, log(-1), log(-1), log(-1), log(-1), log(-1), 0, date);
		
		powspe_new=new double[C_size_powspe];
		dzeros(powspe_new, C_size_powspe);

		for (unsigned int i_gen = 1; i_gen <= C_N_Gen; i_gen++)
		{
/* T */			cout << "Gen: " << i_gen << endl;
			rho_k_gen=new double[C_rho_k_max*C_N_neuron];
			mu_gen=new double[C_N_neuron];
			rate_gen=new double[C_N_neuron];
			sigma_gen=new double[C_N_neuron];
// mu_eps_testaverage(			eps_gen=new double[C_N_neuron];
			for (unsigned int i_neuron = 0; i_neuron < C_N_neuron; i_neuron++)
			{
/* T */				cout << "Neuron: " << i_neuron << endl;			
			/* generate Random-Numbers&Diffusion-Current*/
				rand_id++;
				calc_I_diff(I_input, powspe_old, C_ndt, C_dt, C_rate, C_tau_s, plan_I_input, rand_id, now);

			/* define mu (and eps_diff) + create ISI-train */
				if (arg_mu_cv)
				{
					mu_eps = mu_eps_test(C_rate, C_cv, C_tau_r__dt, C_T_mu, C_tol_mu, C_dt, C_ndt_mu, I_input);
					if (mu_eps[1]==log(-1)) {
						i_neuron--;
						continue;
					}
					interval.lif_neuron(mu_eps[1], C_rate, mu_eps[2], C_dt, C_tau_r__dt, C_ndt, I_input);
				}
				else
				{
					mu = mutest(C_rate, C_eps, C_tau_r__dt, C_T_mu, C_tol_mu, C_dt, C_ndt_mu, I_input);
					interval.lif_neuron(mu, C_rate, C_eps, C_dt, C_tau_r__dt, C_ndt, I_input);
				}

			/* calculate Powerspectrum - normalized to the rate */
				S_temp= new double[C_size_powspe];
				powerspectrum(S_temp, interval, C_ndt, C_dt);

				for (unsigned int i_S = 0; i_S < C_size_powspe; i_S++)
				{
					powspe_new[i_S]+=S_temp[i_S]/C_N_neuron;
				}			
				delete[] S_temp;

			/* calculate averages of rate and mu */
				rate_gen[i_neuron]=interval.rate();
				sigma_gen[i_neuron]=interval.var(C_dt);
				mu_gen[i_neuron]= (arg_mu_cv ? mu_eps[1] : mu);
				interval.calc_rho_k(C_rho_k_max,C_dt, rho_k_gen, i_neuron, interval);
// mu_eps_test				eps_gen [i_neuron]= (arg_mu_cv ? mu_eps[2] : C_eps);
			}

//* T */			time=clock()-tstart;
//* T */			cout << "Time im msec: " << (1000*time/CLOCKS_PER_SEC) << endl;

		/* saving powerspectrum, rate, sigma, mu and constants to file */
			for (unsigned int k = 0; k < C_rho_k_max; k++)
			{
				rho_avg[k]=0;	
				rho_var[k]=0;
				for (unsigned int i = 0; i < C_N_neuron; i++)
				{
					rho_avg[k]+=rho_k_gen[k+i*C_rho_k_max]/C_N_neuron;
				}
				for (unsigned int i = 0; i < C_N_neuron; i++)
				{
					rho_var[k]+=pow((rho_avg[k]-rho_k_gen[k+i*C_rho_k_max]),2)/C_N_neuron;
				}
			}
			rate_avg=average(rate_gen,C_N_neuron);
			sigma_avg=average(sigma_gen,C_N_neuron);
			mu_avg=average(mu_gen,C_N_neuron);
// mu_eps_test			eps_avg=average(eps_gen,C_N_neuron);

			safe_powspe(powspe_new, rho_avg, rho_var, rate_avg, sigma_avg, mu_avg, /*eps_avg,*/variance(rate_gen,C_N_neuron,rate_avg), variance(sigma_gen,C_N_neuron,sigma_avg), variance(mu_gen,C_N_neuron,mu_avg), /*variance(eps_gen,C_N_neuron,eps_avg),*/ i_gen, date);
		
		/* resetting variables */
			delete[] rho_k_gen;
			delete[] rate_gen;
			delete[] sigma_gen;
			delete[] mu_gen;
// mu_eps_test			delete[] eps_gen;
			cpNcl_d_arr(powspe_new, powspe_old, C_size_powspe);
		}
	delete[] I_input;
	delete[] powspe_new;
	delete[] powspe_old;
	fftw_destroy_plan(plan_I_input);
	fftw_forget_wisdom();

	return 0;
}
