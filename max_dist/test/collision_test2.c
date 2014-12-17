/// \file collision_test2.c
/// \brief Test prediction of collisions, at pressure p = 0.
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


static const char *COLL_STR[] = { "COLL_STRING", "COLL_CORE" };


int main()
{
	unsigned int i, im, iw;

	const double pressure = 0;			// pressure is zero
	const double beta = 0.5;			// inverse temperature

	const double mass[2] = { 1, 3 };	// alternating masses

	const unsigned int nsweeps = 13;	// number of 'sweeps' of the event table

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	printf("pressure: %g\n", pressure);
	printf("beta:     %g\n", beta);

	printf("alternating masses: (%g, %g)\n", mass[0], mass[1]);

	printf("N: %d\n", NUM_SITES);
	assert(NUM_SITES == 16);

	// artificial UNIX time
	unsigned int itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// lattice field
	field_t f;
	Field_Create(mass, &f);
	// initialize elongations and momenta using single-site probability density, Eq. (2.11)
	EquilibrateCanonical(pressure, beta, &seed, &f);

	// initial mean values
	double r_mean0 = 0;		// mean of r's (elongations), remains constant in time
	double p_mean0 = 0;		// mean momentum
	double e_mean0 = 0;		// mean energy
	for (im = 0; im < 2; im++)
	{
		for (i = im; i < NUM_SITES; i += 2)
		{
			const lattpoint_t *point = &f.pts[i];
			r_mean0 += point->r;							// elongation
			p_mean0 += mass[im] * point->v;					// momentum
			e_mean0 += LattPoint_Energy(mass[im], point);	// energy
		}
	}
	r_mean0 /= NUM_SITES;
	p_mean0 /= NUM_SITES;
	e_mean0 /= NUM_SITES;
	//e_mean0 -= 0.5 * square(p_mean0);	// internal energy per particle

	// stochastic expectation values
	double r_avr = AverageElongation(pressure, beta);
	double p_avr = 0;
	double e_avr = AverageEnergy(beta);

	// compare mean values with stochastic expectation values
	printf("elongation mean: %g, stochastic average: %g\n", r_mean0, r_avr);
	printf("momentum   mean: %g, stochastic average: %g\n", p_mean0, p_avr);
	printf("energy     mean: %g, stochastic average: %g\n", e_mean0, e_avr);

	// save elongations and momenta to disk
	{
		double r[NUM_SITES];
		double p[NUM_SITES];

		for (im = 0; im < 2; im++)
		{
			for (i = im; i < NUM_SITES; i += 2)
			{
				r[i] = f.pts[i].r;
				p[i] = mass[im] * f.pts[i].v;
			}
		}

		WriteData("../test/collision_test2_field_r0.dat", r, sizeof(double), NUM_SITES, false);
		WriteData("../test/collision_test2_field_p0.dat", p, sizeof(double), NUM_SITES, false);
	}

	// event table
	event_table_t tab;
	assert(NUM_SLOTS == 4);
	EventTable_Allocate(0.1, &tab);	// dt, table
	printf("event table number of slots: %i, slot interval dt: %g, events per slot: %i\n", NUM_SLOTS, tab.dt, EVENTS_PER_SLOT);
	// time period
	const double period = NUM_SLOTS * tab.dt;
	printf("event table time period: %g\n", period);

	// initial kinetic and potential energy
	double en0 = Field_Energy(&f);
	printf("initial energy: %g\n", en0);

	// collect collision events
	CollectEvents(&f, &tab);

	for (iw = 0; iw < nsweeps; iw++)	// several sweeps
	{
		// process events
		SweepTable(&f, &tab);

		// synchronize up to period of event table
		Field_Synchronize(period, &f);

		// subtract 'period' from all timing events
		GlobalTimeshift(period, &f, &tab);
		// correspondingly add 'period' to time offset of event log
		#ifdef _DEBUG
		event_log.time_offset += period;
		#endif
	}

	printf("finished simulation, simulation time: %g\n", nsweeps * period);

	// final mean values (should be conserved)
	double r_mean1 = 0;		// mean of r's (elongations), remains constant in time
	double p_mean1 = 0;		// mean momentum
	double e_mean1 = 0;		// mean energy
	for (im = 0; im < 2; im++)
	{
		for (i = im; i < NUM_SITES; i += 2)
		{
			const lattpoint_t *point = &f.pts[i];
			r_mean1 += point->r;							// elongation
			p_mean1 += mass[im] * point->v;					// momentum
			e_mean1 += LattPoint_Energy(mass[im], point);	// energy
		}
	}
	r_mean1 /= NUM_SITES;
	p_mean1 /= NUM_SITES;
	e_mean1 /= NUM_SITES;
	//e_mean1 -= 0.5 * square(p_mean1);	// internal energy per particle

	// compare initial and final mean values (should be conserved)
	printf("final elongation mean: %g, difference to initial value: %g (should be zero)\n", r_mean1, fabs(r_mean0 - r_mean1));
	printf("final momentum   mean: %g, difference to initial value: %g (should be zero)\n", p_mean1, fabs(p_mean0 - p_mean1));
	printf("final energy     mean: %g, difference to initial value: %g (should be zero)\n", e_mean1, fabs(e_mean0 - e_mean1));

	// final energy
	double en1 = Field_Energy(&f);
	printf("final energy: %g\n", en1);
	printf("energy difference: %g (should be zero)\n", en1 - en0);

	#ifdef _DEBUG

	// report events
	i = 0;
	event_node_t *node = event_log.first;
	double *evtime = malloc(event_log.num * sizeof(double));
	int    *eviprt = malloc(event_log.num * sizeof(int));
	int    *evtype = malloc(event_log.num * sizeof(int));
	while (node != NULL)
	{
		printf("next event: time: %g, iprt: %i, type: %s\n", node->ev.time, node->ev.iprt, COLL_STR[node->ev.type]);
		evtime[i] = node->ev.time;
		eviprt[i] = node->ev.iprt;
		evtype[i] = node->ev.type;
		node = node->next;
		i++;
	}
	assert(i == event_log.num);

	// save event data to disk
	WriteData("../test/collision_test2_event_time.dat", evtime, sizeof(double), event_log.num, false);
	WriteData("../test/collision_test2_event_iprt.dat", eviprt, sizeof(int),    event_log.num, false);
	WriteData("../test/collision_test2_event_type.dat", evtype, sizeof(int),    event_log.num, false);

	// read reference data from disk
	double *evtime_ref = malloc(event_log.num * sizeof(double));
	int    *eviprt_ref = malloc(event_log.num * sizeof(int));
	int    *evtype_ref = malloc(event_log.num * sizeof(int));
	ReadData("../test/collision_test2_ref_event_time.dat", evtime_ref, sizeof(double), event_log.num);
	ReadData("../test/collision_test2_ref_event_iprt.dat", eviprt_ref, sizeof(int),    event_log.num);
	ReadData("../test/collision_test2_ref_event_type.dat", evtype_ref, sizeof(int),    event_log.num);

	// calculate average error
	printf("comparing with reference collision data...\n");
	double err = 0;
	for (i = 0; i < event_log.num; i++)
	{
		err += fabs(evtime[i] - evtime_ref[i]);
		err +=  abs(eviprt[i] - eviprt_ref[i]);
		err +=  abs(evtype[i] - evtype_ref[i]);
	}
	err /= event_log.num;
	printf("average error: %g\n", err);

	// clean up
	free(evtype);
	free(eviprt);
	free(evtime);
	free(evtype_ref);
	free(eviprt_ref);
	free(evtime_ref);
	EventLog_Delete(&event_log);
	#endif	// _DEBUG

	// clean up
	EventTable_Delete(&tab);

	return 0;
}
