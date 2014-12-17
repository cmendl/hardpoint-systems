/// \file time_correlation.c
/// \brief Calculate time correlations using FFT.
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
#include <stdlib.h>
#include <memory.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for time correlation structure
///
void TimeCorrelation_Allocate(const unsigned int num_t, time_correlation_t *tcorr)
{
	int i, j;

	tcorr->num_t = num_t;

	// store averages for twice the time correlation length
	tcorr->avr = fftw_malloc(2*num_t * sizeof(pointdata_t));
	tcorr->sqr = fftw_malloc(  num_t * sizeof(cov_mat_t));

	memset(tcorr->avr, 0, 2*num_t * sizeof(pointdata_t));
	memset(tcorr->sqr, 0,   num_t * sizeof(cov_mat_t));

	// twice the time correlation length to accomodate all time differences
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 3; j++)
		{
			FFTtransform_Allocate(2*num_t, FFTW_FORWARD, &tcorr->fft_trans[i][j]);
		}
	}
	FFTtransform_Allocate(2*num_t, FFTW_BACKWARD, &tcorr->fft_back);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete time correlation structure (free memory)
///
void TimeCorrelation_Delete(time_correlation_t *tcorr)
{
	int i, j;

	FFTtransform_Delete(&tcorr->fft_back);

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 3; j++)
		{
			FFTtransform_Delete(&tcorr->fft_trans[i][j]);
		}
	}

	fftw_free(tcorr->sqr);
	fftw_free(tcorr->avr);

	tcorr->num_t = 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Reset time correlation structure (set all entries to zero)
///
void TimeCorrelation_Reset(time_correlation_t *tcorr)
{
	memset(tcorr->avr, 0, 2*tcorr->num_t * sizeof(pointdata_t));
	memset(tcorr->sqr, 0,   tcorr->num_t * sizeof(cov_mat_t));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Complex multiplication of two numbers 'a' and 'b'
///
static inline void ComplexMultiply(const fftw_complex a, const fftw_complex b, fftw_complex ret)
{
	ret[0] = a[0]*b[0] - a[1]*b[1];
	ret[1] = a[0]*b[1] + a[1]*b[0];
}


//_______________________________________________________________________________________________________________________
///
/// \brief Add 'data' correlations to 'tcorr'
///
/// \param data array of length 2*num_t (twice the time correlation length)
/// \param tcorr correlations are accumulated in 'tcorr'
///
void TimeCorrelation_Collect(const pointdata_t *data, time_correlation_t *tcorr)
{
	unsigned int n;
	int j, k;

	// accumulate averages
	for (n = 0; n < 2*tcorr->num_t; n++)
	{
		tcorr->avr[n].r += data[n].r;
		tcorr->avr[n].p += data[n].p;
		tcorr->avr[n].e += data[n].e;
	}

	// copy data for FFT
	for (j = 0; j < 3; j++)
	{
		// n == 0
		{
			tcorr->fft_trans[0][j].fft_inp[0][0] = data[0].g[j];
			tcorr->fft_trans[1][j].fft_inp[0][0] = data[0].g[j];
		}
		for (n = 1; n < tcorr->num_t; n++)
		{
			// set real parts only
			tcorr->fft_trans[0][j].fft_inp[               n][0] = data[n].g[j];
			tcorr->fft_trans[1][j].fft_inp[2*tcorr->num_t-n][0] = data[n].g[j];		// effectively reverse sign of index
		}
		for (n = tcorr->num_t; n < 2*tcorr->num_t; n++)
		{
			// keep zero entries in 'tcorr->fft_trans[0]'
			tcorr->fft_trans[1][j].fft_inp[2*tcorr->num_t-n][0] = data[n].g[j];		// effectively reverse sign of index
		}
	}

	// perform forward FFTs
	for (k = 0; k < 2; k++)
	{
		for (j = 0; j < 3; j++)
		{
			fftw_execute(tcorr->fft_trans[k][j].plan);
		}
	}

	// normalization factor; '2*tcorr->num_t' from FFT
	const double invN2 = 1.0 / (2*tcorr->num_t * tcorr->num_t);

	for (j = 0; j < 3; j++)
	{
		for (k = 0; k < 3; k++)
		{
			// pointwise multiplication
			for (n = 0; n < 2*tcorr->num_t; n++)
			{
				ComplexMultiply(tcorr->fft_trans[0][j].fft_out[n], tcorr->fft_trans[1][k].fft_out[n], tcorr->fft_back.fft_inp[n]);
			}

			// perform backward FFT
			fftw_execute(tcorr->fft_back.plan);

			// accumulate normalized data
			// n == 0
			{
				tcorr->sqr[0].mat[j + 3*k] += tcorr->fft_back.fft_out[0][0] * invN2;
			}
			for (n = 1; n < tcorr->num_t; n++)
			{
				tcorr->sqr[n].mat[j + 3*k] += tcorr->fft_back.fft_out[2*tcorr->num_t - n][0] * invN2;
			}
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Multiply by scalar 'x', i.e., tcorr.avr *= x and tcorr.sqr *= x
///
void TimeCorrelation_ScalarMultiply(const double x, time_correlation_t *tcorr)
{
	unsigned int n;
	int j;

	for (n = 0; n < 2*tcorr->num_t; n++)
	{
		tcorr->avr[n].r *= x;
		tcorr->avr[n].p *= x;
		tcorr->avr[n].e *= x;
	}

	for (n = 0; n < tcorr->num_t; n++)
	{
		for (j = 0; j < 9; j++)
		{
			tcorr->sqr[n].mat[j] *= x;
		}
	}
}
