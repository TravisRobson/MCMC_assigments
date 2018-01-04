#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include "Constants.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage)
{
    double val = (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3.1f%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

void make_data(double *data, long N, gsl_rng *r);
double get_logL_sampler(double x);

int main(int argc, char *argv[])
{
	int seed;
	
	long i;
	const long NDATA = 1000;
//	const long NMCMC = 1000;
//	long Naccept = 0;
	
	double time_spent, accept;
// 	double mu_x, nu_x, sigma_x;
// 	double mu_y, nu_y, sigma_y;

	double *data;
	double *evals, **evecs;
	
	clock_t begin = clock();
	clock_t end;
	
	FILE *out_file;
	
	fprintf(stdout, "==============================================================\n\n");
	
	// set up the GSL random number stuff
	seed = atoi(argv[1]);
	gsl_rng_env_setup();
	
	const gsl_rng_type *TT = gsl_rng_default;
	gsl_rng *r 			   = gsl_rng_alloc(TT);
	gsl_rng_set(r, seed);
	
	out_file = fopen(argv[2], "w");
	
	// construct data
	data = malloc(NDATA*sizeof(double));
	make_data(data, NDATA, r);
	// FILE *dat_file; dat_file = fopen("samples.dat", "w");
	// for (i=0; i<NDATA; i++) fprintf(, "%.12g\n", data[i]); // prints data
	// fclose(dat_file);
	
	// store eigen-bs
	
	
	
	
	
	
	fclose(out_file);
// 	accept = (double)Naccept/(double)NMCMC;
// 	fprintf(stdout, "Acceptance Rate: %.2f%%\n", accept*100);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nMCMC duration: %f sec\n", time_spent);
	
	
	
	fprintf(stdout, "\n==============================================================\n");
	
	
	free(data);
	
	free(evals);
	free(evecs[0]); free(evecs[1]); free(evecs[2]);
	
	return 0;
}

void make_data(double *data, long N, gsl_rng *r)
{		
	const int skip = 100;	// roughly the auto-correlation length for the Fisher steps

	long i, j;
	
	double Fish;
	double x, y;
	double logLy, logLx, loga;
	
	// set the proposal distribution variance
	Fish = 1./sqrt((1. + nu)/(nu*sigma*sigma));
	
	// set the initial samples
	x = mu;
	y = mu;
	
	logLx = get_logL_sampler(x);
	
	for (i=0; i<N; i++)
	{
		// wait for an independent samples
		for (j=0; j<skip; j++)
		{
			y = x + gsl_ran_gaussian(r, Fish);
			logLy = get_logL_sampler(y);
			loga = log(gsl_rng_uniform(r));
		
			if (logLy > -INFINITY && loga < (logLy - logLx))
			{
				x = y;
				logLx = logLy;
			}
		}
		
		data[i] = x;
	}

	return;
}

double get_logL_sampler(double x)
{
	double logL;
	
    logL  = -0.5*(log(M_PI) + log(nu)) - log(sigma) - log(gsl_sf_gamma(0.5*nu));
    logL += log(gsl_sf_gamma(0.5*(1. + nu))) - 0.5*(1. + nu)*log(1. + (x - mu)*(x - mu)/nu/sigma/sigma);
	
	return logL;
}












