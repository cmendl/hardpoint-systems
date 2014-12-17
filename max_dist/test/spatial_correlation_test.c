/// \file spatial_correlation_test.c
/// \brief Test the calculation of spatial (lattice site) correlations.
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
#include "spatial_correlation_ref.h"
#include "equilibration.h"
#include "field.h"
#include "random.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


int main()
{
	unsigned int i, j, k;

	const double pressure = 2.0;		// pressure
	const double beta = 0.5;			// inverse temperature

	const double mass[2] = { 1, 3 };	// alternating masses

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	printf("pressure: %g\n", pressure);
	printf("beta:     %g\n", beta);

	printf("alternating masses: (%g, %g)\n", mass[0], mass[1]);

	printf("N: %d\n", NUM_SITES);
	assert(NUM_SITES == 4096);

	// artificial UNIX time
	unsigned int itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// lattice point data
	pointdata_t *data_latt[2];

	// create two random fields
	for (j = 0; j < 2; j++)
	{
		data_latt[j] = fftw_malloc(NUM_SITES * sizeof(pointdata_t));

		// lattice field
		field_t f;
		Field_Create(mass, &f);

		// initialize elongations and momenta using single-site probability density, Eq. (2.11)
		EquilibrateCanonical(pressure, beta, &seed, &f);

		// store elongation and momentum, and calculate energies of field variables
		for (k = 0; k < 2; k++)
		{
			for (i = k; i < NUM_SITES; i += 2)
			{
				const lattpoint_t *point = &f.pts[i];
				data_latt[j][i].r = point->r;							// elongation
				data_latt[j][i].p = mass[k] * point->v;					// momentum
				data_latt[j][i].e = LattPoint_Energy(mass[k], point);	// energy
			}
		}
	}

	// spatial correlation structure
	spatial_correlation_t scorr;
	SpatialCorrelation_Allocate(&scorr);
	SpatialCorrelation_Collect(data_latt[0], data_latt[1], &scorr);

	// reference spatial correlation structure
	spatial_correlation_ref_t scorr_ref;
	SpatialCorrelationRef_Allocate(&scorr_ref);
	SpatialCorrelationRef_Collect(data_latt[0], data_latt[1], &scorr_ref);

	double err = 0;
	for (i = 0; i < NUM_SITES; i++)
	{
		for (j = 0; j < 9; j++)
		{
			err = fmax(err, fabs(scorr.sqr[i].mat[j] - scorr_ref.sqr[i].mat[j]));
		}
	}
	for (j = 0; j < 2; j++)
	{
		err = fmax(err, fabs(scorr.avr[j].r - scorr_ref.avr[j].r));
		err = fmax(err, fabs(scorr.avr[j].p - scorr_ref.avr[j].p));
		err = fmax(err, fabs(scorr.avr[j].e - scorr_ref.avr[j].e));
	}
	printf("maximum error: %g\n", err);

	// clean up
	SpatialCorrelationRef_Delete(&scorr_ref);
	SpatialCorrelation_Delete(&scorr);
	for (j = 0; j < 2; j++)
	{
		fftw_free(data_latt[j]);
	}
	fftw_cleanup();

	return 0;
}
