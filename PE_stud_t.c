#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include "Constants.h"
#include "IO_ass2.h"


void make_data(double *data, long N, gsl_rng *r);
double get_logL_sampler(double x);
double get_logL(double *data, long N, double *params);
void check_priors(double *params, int *meet_priors);

int main(int argc, char *argv[])
{
	int seed;
	const int NP = 3;
	int j, k;
	int meet_priors = 1;
	
	long i;
	
	const long NDATA = 1000;
	const long NMCMC = 1e5;
	long Naccept = 0;
	
	double time_spent, accept, jump;
	double logLx, logLy, loga, logL_max;
 	double *params_x, *params_y, *params_max;

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
	evals = malloc(NP*sizeof(double));
	evecs = malloc(NP*sizeof(double *));
	for (i=0; i<NP; i++) evecs[i] = malloc(NP*sizeof(double));
	
	evals[0] = (7. + 3.*pow(M_PI,2.) + sqrt(4261. - 390.*pow(M_PI,2.) + 9.*pow(M_PI,4.)))/72.;
	evals[1] = 0.6666666666666666;
	evals[2] = (7. + 3.*pow(M_PI,2.) - sqrt(4261. - 390.*pow(M_PI,2.) + 9.*pow(M_PI,4.)))/72.;
	
	evecs[0][0] = 0.;
	evecs[0][1] = (65. - 3.*pow(M_PI,2.) - sqrt(4261. - 390.*pow(M_PI,2.) + 9.*pow(M_PI,4.)))/6.;
	evecs[0][2] = 1.;
	
	evecs[1][0] = 1.;
	evecs[1][1] = 0.;
	evecs[1][2] = 0.;
	
	evecs[2][0] = 0.;
	evecs[2][1] = (65. - 3.*pow(M_PI,2.) + sqrt(4261. - 390.*pow(M_PI,2) + 9.*pow(M_PI,4.)))/6.;
	evecs[2][2] = 1.;
	
	params_x   = malloc(NP*sizeof(double));
	params_y   = malloc(NP*sizeof(double));
	params_max = malloc(NP*sizeof(double));
	
	for (i=0; i<NP; i++)
	{
		params_x[i]   = 0.;
		params_y[i]   = 0.;
		params_max[i] = 0.;
	}
	
	// start at the true values
	params_x[0] = mu;
	params_x[1] = nu;
	params_x[2] = sigma;
	
	logLx = get_logL(data, NDATA, params_x);
	logLy = logLx;
	logL_max = -1.0e30;
	
	fprintf(stdout, "initial logL: %f\n", logLx);
	
	for (i=0; i<NMCMC; i++)
	{
		// propose new parameters
		k = (int)(NP*gsl_rng_uniform(r));
		for (j=0; j<NP; j++)
		{
			jump = 0.1*gsl_ran_gaussian(r, 1.)*evecs[k][j]/sqrt(evals[k]*NP);
			params_y[j] = params_x[j] + jump;
		}
		check_priors(params_y, &meet_priors);

		if (meet_priors == 1) logLy = get_logL(data, NDATA, params_y);
		loga = log(gsl_rng_uniform(r));
		
		if (logLy > -INFINITY && loga < (logLy - logLx) && meet_priors == 1)
		{
			Naccept++;
			
			for (j=0; j<NP; j++) params_x[j] = params_y[j];
			logLx = logLy;
			
			// check for max in posterior
			if (logLx > logL_max)
			{	
				logL_max = logLx;
				for (j=0; j<NP; j++) params_max[j] = params_x[j];
			}
		}
		
 		fprintf(out_file, "%.12g %.12g %.12g %.12g\n", logLx, params_x[0], params_x[1], params_x[2]);
		meet_priors = 1; //reset
	}
	
	
	
	fclose(out_file);
	accept = (double)Naccept/(double)NMCMC;
	fprintf(stdout, "Acceptance Rate: %.2f%%\n", accept*100);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nMCMC duration: %f sec\n\n", time_spent);
	
	fprintf(stdout, "logL_max: %f\n", logL_max);
	fprintf(stdout, "-------------------------------\n");
	fprintf(stdout, "max mu: %f\n", params_max[0]);
	fprintf(stdout, "max nu: %f\n", params_max[1]);
	fprintf(stdout, "max sg: %f\n", params_max[2]);
	
	
	
	
	fprintf(stdout, "\n==============================================================\n");
	
	
	free(data);
	
	free(params_x);
	free(params_y);
	free(params_max);
	
	free(evals);
	free(evecs[0]); free(evecs[1]); free(evecs[2]);
	
	return 0;
}

void check_priors(double *params, int *meet_priors)
{
	if (params[0] < -5. || params[0] >  5.) *meet_priors = 0;
	if (params[1] < 0.1 || params[1] > 10.) *meet_priors = 0;
	if (params[2] < 0.1 || params[2] > 10.) *meet_priors = 0;
	
	return;
}

double get_logL(double *data, long N, double *params)
{
	long i;
	
	double logL = 0.;
	double a;
	double nu_cur, mu_cur, sg_cur;
	
	mu_cur = params[0];
	nu_cur = params[1];
	sg_cur = params[2];
	
	
	a = -0.5*(log(M_PI) + log(nu_cur)) - log(sg_cur) - log(gsl_sf_gamma(0.5*nu_cur)) 
					    + log(gsl_sf_gamma(0.5*(1. + nu_cur)));
	
	for (i=0; i<N; i++)
	{
		logL += -0.5*(1. + nu_cur)*log(1. + (data[i] - mu_cur)*(data[i] - mu_cur)/nu_cur/sg_cur/sg_cur);
	}
	
	logL += N*a;
	
	return logL;
//	return 1.;
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












