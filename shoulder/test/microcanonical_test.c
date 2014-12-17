/// \file microcanonical_test.c
/// \brief Test microcanonical equilibration.
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

#include "equilibration.h"
#include "spatial_correlation.h"
#include "util.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif



int main()
{
	unsigned int is, ip;

	const double pressure = 1.2;		// pressure
	const double beta = 2.0;			// inverse temperature

	const double length = AverageElongation(pressure, beta);
	const double energy = AverageEnergy(pressure, beta);

	printf("pressure: %g\n", pressure);
	printf("beta:     %g\n", beta);

	printf("average length: %g\n", length);
	printf("average energy: %g\n", energy);

	printf("N: %i\n", NUM_SITES);
	assert(NUM_SITES == 16);

	unsigned int nsamples = 100000;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// artificial UNIX time
	unsigned int itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// lattice field
	field_t f;
	Field_Create(&f);

	// spatial correlation structure
	spatial_correlation_t scorr;
	SpatialCorrelation_Allocate(&scorr);
	// lattice point data
	pointdata_t pointdata[NUM_SITES];

	// start timer
	clock_t t_start = clock();

	for (is = 0; is < nsamples; is++)
	{
		EquilibrateMicrocanonical(length, energy, &seed, &f);

		// store elongation and momentum, and calculate energies of field variables
		for (ip = 0; ip < NUM_SITES; ip++)
		{
			const lattpoint_t *point = &f.pts[ip];
			pointdata[ip].r = point->r;						// elongation
			pointdata[ip].p = point->p;						// momentum
			pointdata[ip].e = LattPoint_Energy(point);		// energy
		}

		// collect and accumulate spatial auto-correlations
		SpatialCorrelation_Collect(pointdata, pointdata, &scorr);
	}

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("finished microcanonical equilibration, CPU time: %g, average CPU time per run: %g\n", cpu_time, cpu_time / (nsamples));

	// divide by 'nsamples'
	SpatialCorrelation_ScalarMultiply(1.0 / nsamples, &scorr);

	// save correlations to disk
	WriteData("../test/microcanonical_test_avr0.dat", &scorr.avr[0], sizeof(pointdata_t), 1, false);
	WriteData("../test/microcanonical_test_avr1.dat", &scorr.avr[1], sizeof(pointdata_t), 1, false);
	WriteData("../test/microcanonical_test_sqrs.dat", scorr.sqr, sizeof(cov_mat_t), NUM_SITES, false);

	// clean up
	SpatialCorrelation_Delete(&scorr);
	fftw_cleanup();

	return 0;
}
