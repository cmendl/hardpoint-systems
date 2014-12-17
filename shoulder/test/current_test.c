/// \file current_test.c
/// \brief Test current calculation.
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

#include "simulation.h"
#include "equilibration.h"
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
	unsigned int i;

	const double pressure = 1.2;		// pressure
	const double beta = 2.0;			// inverse temperature

	const unsigned int nsweeps = 5;		// number of 'sweeps' of the event table

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	printf("pressure: %g\n", pressure);
	printf("beta:     %g\n", beta);

	printf("N: %d\n", NUM_SITES);
	assert(NUM_SITES == 16);

	// canonical ensemble parameters
	ensemble_params_t params;
	EnsembleParams_Fill(pressure, beta, &params);

	// artificial UNIX time
	unsigned int itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// lattice field
	field_t f;
	Field_Create(&f);
	// initialize elongations and momenta using single-site probability density, Eq. (2.11)
	EquilibrateCanonical(&params, &seed, &f);

	// save elongations and momenta to disk
	{
		double r[NUM_SITES];
		double p[NUM_SITES];
		for (i = 0; i < NUM_SITES; i++)
		{
			r[i] = f.pts[i].r;
			p[i] = f.pts[i].p;
		}

		WriteData("../test/current_test_field_r0.dat", r, sizeof(double), NUM_SITES, false);
		WriteData("../test/current_test_field_p0.dat", p, sizeof(double), NUM_SITES, false);
	}

	// event table
	event_table_t tab;
	assert(NUM_SLOTS == 4);
	EventTable_Allocate(0.1, &tab);		// dt, table
	// time period
	const double period = NUM_SLOTS * tab.dt;
	printf("event table time period : %g\n", period);

	// momentum and energy current for each sweep
	current_t *J = malloc(nsweeps * sizeof(current_t));

	// initial kinetic and potential energy
	double en0 = Field_Energy(&f);
	printf("initial energy: %g\n", en0);

	// collect collision events
	CollectEvents(&f, &tab);

	for (i = 0; i < nsweeps; i++)	// several sweeps
	{
		// process events and record current
		J[i] = SweepTable(&f, &tab);

		printf("momentum and energy current during sweep %i: (%g, %g)\n", i, J[i].p, J[i].e);

		// synchronize up to period of event table
		//Field_Synchronize(period, &f);

		// subtract 'period' from all timing events
		GlobalTimeshift(period, &f, &tab);
		// correspondingly add 'period' to time offset of event log
		#ifdef _DEBUG
		event_log.time_offset += period;
		#endif
	}

	printf("finished simulation, simulation time: %g\n", nsweeps * period);

	// final kinetic and potential energy
	double en1 = Field_Energy(&f);
	printf("final energy: %g\n", en1);
	printf("energy difference: %g (should be zero)\n", en1 - en0);

	// read reference current data from disk
	current_t *J_ref = malloc(nsweeps * sizeof(current_t));
	ReadData("../test/current_test_Jref.dat", J_ref, sizeof(current_t), nsweeps);

	// calculate average error
	printf("comparing with reference data...\n");
	current_t err = { 0 };
	for (i = 0; i < nsweeps; i++)
	{
		err.p += fabs(J[i].p - J_ref[i].p);
		err.e += fabs(J[i].e - J_ref[i].e);
	}
	err.p /= nsweeps;
	err.e /= nsweeps;
	printf("average momentum and energy current error: (%g, %g)\n", err.p, err.e);

	#ifdef _DEBUG
	EventLog_Delete(&event_log);
	#endif	// _DEBUG

	// clean up
	free(J_ref);
	free(J);
	EventTable_Delete(&tab);

	return 0;
}
