/// \file spatial_correlation_ref.c
/// \brief Calculate spatial (lattice site) correlations, reference implementation for testing.
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

#include "spatial_correlation_ref.h"
#include "field.h"
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for spatial correlation structure
///
void SpatialCorrelationRef_Allocate(spatial_correlation_ref_t *scorr)
{
	//scorr->N = N;
	memset(scorr->avr, 0, sizeof(scorr->avr));
	scorr->sqr = calloc(NUM_SITES, sizeof(cov_mat_t));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete spatial correlation structure (free memory)
///
void SpatialCorrelationRef_Delete(spatial_correlation_ref_t *scorr)
{
	free(scorr->sqr);
	//scorr->N = 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Reset spatial correlation structure (set all entries to zero)
///
void SpatialCorrelationRef_Reset(spatial_correlation_ref_t *scorr)
{
	memset(scorr->avr, 0, sizeof(scorr->avr));
	memset(scorr->sqr, 0, NUM_SITES * sizeof(cov_mat_t));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Outer product of lattice point data
///
static inline void Data_Product(const pointdata_t *avr0, const pointdata_t *avr1, cov_mat_t *prod)
{
	prod->entries.rr = avr0->r * avr1->r;
	prod->entries.pr = avr0->p * avr1->r;
	prod->entries.er = avr0->e * avr1->r;

	prod->entries.rp = avr0->r * avr1->p;
	prod->entries.pp = avr0->p * avr1->p;
	prod->entries.ep = avr0->e * avr1->p;

	prod->entries.re = avr0->r * avr1->e;
	prod->entries.pe = avr0->p * avr1->e;
	prod->entries.ee = avr0->e * avr1->e;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Fill spatial correlation structure
///
void SpatialCorrelationRef_Collect(const pointdata_t *data0, const pointdata_t *data1, spatial_correlation_ref_t *scorr)
{
	int i, j, k;

	// initially set to zero
	memset(scorr->sqr, 0, NUM_SITES * sizeof(cov_mat_t));
	memset(scorr->avr, 0, sizeof(scorr->avr));

	for (i = 0; i < NUM_SITES; i++)
	{
		// reference lattice point data
		const pointdata_t *dt0 = &data0[i];

		for (j = 0; j < NUM_SITES; j++)	// offset delta
		{
			cov_mat_t prod;
			Data_Product(dt0, &data1[(i + j) & SITE_MASK], &prod);

			// accumulate
			for (k = 0; k < 9; k++) {
				scorr->sqr[j].mat[k] += prod.mat[k];
			}
		}

		scorr->avr[0].r += dt0->r;
		scorr->avr[0].p += dt0->p;
		scorr->avr[0].e += dt0->e;

		// reference lattice point data
		const pointdata_t *dt1 = &data1[i];

		scorr->avr[1].r += dt1->r;
		scorr->avr[1].p += dt1->p;
		scorr->avr[1].e += dt1->e;
	}

	// divide by N
	const double invN = 1.0 / NUM_SITES;

	scorr->avr[0].r *= invN;
	scorr->avr[0].p *= invN;
	scorr->avr[0].e *= invN;

	scorr->avr[1].r *= invN;
	scorr->avr[1].p *= invN;
	scorr->avr[1].e *= invN;

	for (j = 0; j < NUM_SITES; j++)
	{
		for (k = 0; k < 9; k++)
		{
			scorr->sqr[j].mat[k] *= invN;
		}
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Multiply by scalar 'x', i.e., scorr *= x
///
void SpatialCorrelationRef_ScalarMultiply(const double x, spatial_correlation_ref_t *scorr)
{
	int i, j;

	scorr->avr[0].r *= x;
	scorr->avr[0].p *= x;
	scorr->avr[0].e *= x;

	scorr->avr[1].r *= x;
	scorr->avr[1].p *= x;
	scorr->avr[1].e *= x;

	for (i = 0; i < NUM_SITES; i++) {
		for (j = 0; j < 9; j++) {
			scorr->sqr[i].mat[j] *= x;
		}
	}
}
