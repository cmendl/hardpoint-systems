/// \file spatial_correlation.c
/// \brief Calculate spatial (lattice site) correlations using FFT.
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

#include "spatial_correlation.h"
#include "field.h"
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for spatial correlation structure
///
void SpatialCorrelation_Allocate(spatial_correlation_t *scorr)
{
	int i, j;

	//scorr->N = N;
	scorr->sqr = fftw_malloc(NUM_SITES * sizeof(cov_mat_t));
	memset(scorr->avr, 0, sizeof(scorr->avr));
	memset(scorr->sqr, 0, NUM_SITES * sizeof(cov_mat_t));

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 3; j++)
		{
			FFTtransform_Allocate(NUM_SITES, FFTW_FORWARD, &scorr->fft_trans[i][j]);
		}
	}
	FFTtransform_Allocate(NUM_SITES, FFTW_BACKWARD, &scorr->fft_back);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete spatial correlation structure (free memory)
///
void SpatialCorrelation_Delete(spatial_correlation_t *scorr)
{
	int i, j;

	FFTtransform_Delete(&scorr->fft_back);

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 3; j++)
		{
			FFTtransform_Delete(&scorr->fft_trans[i][j]);
		}
	}

	fftw_free(scorr->sqr);

	//scorr->N = 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Reset spatial correlation structure (set all entries to zero)
///
void SpatialCorrelation_Reset(spatial_correlation_t *scorr)
{
	memset(scorr->avr, 0, sizeof(scorr->avr));
	memset(scorr->sqr, 0, NUM_SITES * sizeof(cov_mat_t));
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
/// \brief Fill spatial correlation structure
///
void SpatialCorrelation_Collect(const pointdata_t *data0, const pointdata_t *data1, spatial_correlation_t *scorr)
{
	unsigned int i;
	int j, k;

	// mean values
	pointdata_t avr[2] = { 0 };
	for (i = 0; i < NUM_SITES; i++)
	{
		// reference lattice point data
		const pointdata_t *dt0 = &data0[i];
		const pointdata_t *dt1 = &data1[i];

		for (j = 0; j < 3; j++)
		{
			avr[0].g[j] += dt0->g[j];
			avr[1].g[j] += dt1->g[j];
		}
	}
	// divide by N
	const double invN = 1.0 / NUM_SITES;
	for (j = 0; j < 3; j++)
	{
		avr[0].g[j] *= invN;
		avr[1].g[j] *= invN;
	}

	// accumulate averages
	for (j = 0; j < 3; j++)
	{
		scorr->avr[0].g[j] += avr[0].g[j];
		scorr->avr[1].g[j] += avr[1].g[j];
	}

	// copy data for FFT
	for (i = 0; i < NUM_SITES; i++)
	{
		// reference lattice point data
		const pointdata_t *dt0 = &data0[i];
		const pointdata_t *dt1 = &data1[(NUM_SITES - i) & SITE_MASK];	// effectively reverse sign of index

		// set real parts only
		for (j = 0; j < 3; j++)
		{
			scorr->fft_trans[0][j].fft_inp[i][0] = dt0->g[j];
			scorr->fft_trans[1][j].fft_inp[i][0] = dt1->g[j];
		}
	}

	// perform forward FFTs
	for (k = 0; k < 2; k++)
	{
		for (j = 0; j < 3; j++)
		{
			fftw_execute(scorr->fft_trans[k][j].plan);
		}
	}

	// divide by N^2 (one factor of N from FFT)
	const double invN2 = 1.0 / (NUM_SITES * NUM_SITES);

	for (j = 0; j < 3; j++)
	{
		for (k = 0; k < 3; k++)
		{
			// pointwise multiplication
			for (i = 0; i < NUM_SITES; i++)
			{
				ComplexMultiply(scorr->fft_trans[0][j].fft_out[i], scorr->fft_trans[1][k].fft_out[i], scorr->fft_back.fft_inp[i]);
			}

			// perform backward FFT
			fftw_execute(scorr->fft_back.plan);

			// accumulate normalized data
			for (i = 0; i < NUM_SITES; i++)
			{
				scorr->sqr[i].mat[j + 3*k] += scorr->fft_back.fft_out[(NUM_SITES - i) & SITE_MASK][0] * invN2;
			}
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Multiply by scalar 'x', i.e., scorr.avr *= x and scorr.sqr *= x
///
void SpatialCorrelation_ScalarMultiply(const double x, spatial_correlation_t *scorr)
{
	unsigned int i;
	int j;

	for (j = 0; j < 3; j++)
	{
		scorr->avr[0].g[j] *= x;
		scorr->avr[1].g[j] *= x;
	}

	for (i = 0; i < NUM_SITES; i++)
	{
		for (j = 0; j < 9; j++)
		{
			scorr->sqr[i].mat[j] *= x;
		}
	}
}
