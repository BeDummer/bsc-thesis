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
void whitenoise_PS(double *arr, const int n, const double r_0) // generates bandlimited white-noise powerspectrum with rate r_0 und cut-off frequency f_c=1/(2*C_dt)
{
	for (unsigned int i=0;i<n;i++)
		arr[i]=r_0;
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
int main()
{
/* variables for calculation */
	double powspe_old[C_size_powspe]={0.0};
	double powspe_new[C_size_powspe]={0.0};
	double* S_temp;
	double* I_input;
	double mu, mu_gen=0., rate_gen=0.;
	ISI interval(C_T_max, C_rate);

/* variables for initialization of random number generators */	
	unsigned long int rand_id=0;
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
		whitenoise_PS(powspe_old,C_size_powspe,C_rate);
		dzeros(powspe_new, C_size_powspe);

		for (unsigned int i_gen = 0; i_gen < C_N_Gen; i_gen++)
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
				mu_gen+=mu/C_N_neuron;
			}

//* T */			time=clock()-tstart;
//* T */			cout << "Time im msec: " << (1000*time/CLOCKS_PER_SEC) << endl;

		/* write Powerspectrum, mu and rate into buffer */
			buffer << "dt\t" << "N\t" << "eps\t" << "tau_r\t" << "N_neuron\n" << C_dt << "\t" << C_ndt << "\t" << "\t" << C_eps <<"\t" << C_tau_r << "\t" << C_N_neuron << "\n\nmu\t\trate\n" << mu_gen << "\t" << rate_gen << "\n\n";

			for (unsigned int i_safe = 0; i_safe < C_size_powspe; i_safe++)
			{
				buffer << (powspe_new[i_safe]/rate_gen) << endl;
			}

			filename << "data/" << date << "__" << i_gen << ".dat";
			filename_tmp=filename.str();

		/* saving data to file from buffer */
			file.open(filename_tmp.c_str());
				file << buffer.str();
			file.close();

		/* resetting variables */
			filename.str("");
			buffer.str("");
			cpNcl_d_arr(powspe_new, powspe_old, C_size_powspe);
			rate_gen=0.;
			mu_gen=0.;
		}
	delete[] I_input;
	fftw_destroy_plan(plan_I_input);
	fftw_forget_wisdom();

	return 0;
}
