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

double get_logL(double x);

int main(int argc, char *argv[])
{
	int seed;
	
	long i;
	const long NMCMC = 1000;
	long Naccept = 0;
	
	double time_spent;
	double alpha;
	double x, y;
	double logLx, logLy, loga;
	double accept;
	double Fish;
	
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
	
	out_file = fopen(argv[3], "w");
	
	// set the proposal distribution variance
	alpha = atof(argv[2]);
	Fish = 1./sqrt((1. + nu)/(nu*sigma*sigma));
	
	// set the initial samples
	x = mu;
	y = mu;
	
	logLx = get_logL(x);
	
	for (i=0; i<NMCMC; i++)
	{
		//y = x + gsl_ran_gaussian(r, alpha);
		y = x + gsl_ran_gaussian(r, Fish);
		logLy = get_logL(y);
		loga = log(gsl_rng_uniform(r));
		
		if (logLy > -INFINITY && loga < (logLy - logLx))
		{
			Naccept++;
			
			x = y;
			logLx = logLy;
		}
		
		fprintf(out_file, "%.12g %.12g\n", logLx, x);
	}
	fclose(out_file);
	
	accept = (double)Naccept/(double)NMCMC;
	fprintf(stdout, "Acceptance Rate: %.2f%%\n", accept*100);
	
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	fprintf(stdout, "\nMCMC duration: %f sec\n", time_spent);
	
	
	
	fprintf(stdout, "\n==============================================================\n");
	
	
	return 0;
}

double get_logL(double x)
{
	double logL;
	
    logL  = -0.5*(log(M_PI) + log(nu)) - log(sigma) - log(gsl_sf_gamma(0.5*nu));
    logL += log(gsl_sf_gamma(0.5*(1. + nu))) - 0.5*(1. + nu)*log(1. + (x - mu)*(x - mu)/nu/sigma/sigma);
	
	return logL;
}












