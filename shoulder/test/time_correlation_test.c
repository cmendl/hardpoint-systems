/// \file time_correlation_test.c
/// \brief Test the calculation of time correlations.
//
//	Copyright (c) 2014, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	References:
//	- Christian B. Mendl, Herbert Spohn
//	  Current fluctuations for anharmonic chains in thermal equilibrium
//	  arXiv:1412.4609
//
//	- Christian B. Mendl, Herbert Spohn
//	  Equilibrium time-correlation functions for one-dimensional hard-point systems
//	  Phys. Rev. E 90, 012147 (2014), arXiv:1403.0213
//
//	- Christian B. Mendl, Herbert Spohn
//	  Dynamic correlators of Fermi-Pasta-Ulam chains and nonlinear fluctuating hydrodynamics
//	  Phys. Rev. Lett. 111, 230601 (2013), arXiv:1305.1209
//
//	- Herbert Spohn
//	  Nonlinear fluctuating hydrodynamics for anharmonic chains
//	  J. Stat. Phys. 154, 1191-1227 (2014), arXiv:1305.6412
//_______________________________________________________________________________________________________________________
//

#include "time_correlation.h"
#include "util.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


// reference implementation
void CalculateAveragedTimeCorrelation(const unsigned int num_t, const pointdata_t *data, cov_mat_t *corr)
{
	unsigned int n, m;
	int j, k;

	for (k = 0; k < 3; k++)
	{
		for (j = 0; j < 3; j++)
		{
			for (m = 0; m < num_t; m++)
			{
				double c = 0;
				for (n = 0; n < num_t; n++)
				{
					c += data[n].g[j] * data[n + m].g[k];
				}
				// normalization
				c /= num_t;

				corr[m].mat[j + 3*k] = c;
			}
		}
	}
}


//_______________________________________________________________________________________________________________________
//


int main()
{
	unsigned int n;
	int l, j;
	#define nsmpl 5
	#define num_t 17

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// random generator seed
	unsigned int iseed = (unsigned int)time(NULL);
	srand(iseed);

	pointdata_t data[nsmpl][2*num_t];
	for (l = 0; l < nsmpl; l++)
	{
		for (n = 0; n < 2*num_t; n++)
		{
			for (j = 0; j < 3; j++)
			{
				data[l][n].g[j] = (2.0 * rand() / RAND_MAX) - 1;
			}
		}
	}

	// calculate reference values
	pointdata_t ref_avr[num_t] = { 0 };
	cov_mat_t   ref_sqr[num_t] = { 0 };
	for (l = 0; l < nsmpl; l++)
	{
		cov_mat_t cur_sqr[num_t];
		CalculateAveragedTimeCorrelation(num_t, data[l], cur_sqr);
		for (n = 0; n < num_t; n++)
		{
			for (j = 0; j < 3; j++)
			{
				ref_avr[n].g[j] += data[l][n].g[j] / nsmpl;
			}
			for (j = 0; j < 9; j++)
			{
				ref_sqr[n].mat[j] += cur_sqr[n].mat[j] / nsmpl;
			}
		}
	}

	time_correlation_t tcorr;
	TimeCorrelation_Allocate(num_t, &tcorr);

	// collect samples
	for (l = 0; l < nsmpl; l++)
	{
		TimeCorrelation_Collect(data[l], &tcorr);
	}
	TimeCorrelation_ScalarMultiply(1.0/nsmpl, &tcorr);

	// compare
	double err_avr = 0;
	double err_sqr = 0;
	for (n = 0; n < num_t; n++)
	{
		for (j = 0; j < 3; j++)
		{
			err_avr = maxf(err_avr, fabs(tcorr.avr[n].g[j] - ref_avr[n].g[j]));
		}
		for (j = 0; j < 9; j++)
		{
			err_sqr = maxf(err_sqr, fabs(tcorr.sqr[n].mat[j] - ref_sqr[n].mat[j]));
		}
	}
	printf("maximum error of averages: %g\n", err_avr);
	printf("maximum error of squares:  %g\n", err_sqr);

	// clean up
	TimeCorrelation_Delete(&tcorr);
	fftw_cleanup();

	return 0;
}
