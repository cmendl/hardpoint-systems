/// \file collision_test.c
/// \brief Test prediction of collisions.
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
#include "random.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


static const char *COLL_STR[] = { "COLL_PLAT_IN", "COLL_PLAT_OUT", "COLL_PLAT_BOUNCE", "COLL_CORE" };


int main()
{
	unsigned int i;
	const double L = 12;		// total length (sum of r's, remains constant in time)
	const double beta = 2.0;	// inverse temperature

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	printf("L:    %g\n", L);
	printf("beta: %g\n", beta);

	printf("N: %d\n", NUM_SITES);
	assert(NUM_SITES == 8);

	// artificial UNIX time
	unsigned int itime = 1420000241;

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// lattice field
	field_t f;
	Field_Create(&f);
	// initialize elongations and momenta (should actually use Eq. (2.11))
	{
		double sigma = 1 / sqrt(beta);
		double ell = L / NUM_SITES;

		for (i = 0; i < NUM_SITES; i++)
		{
			f.pts[i].r = Random_GetUniform(&seed) * ell;
			f.pts[i].p = sigma * Random_GetNormal(&seed);
		}
		// ensure that extensive constraints are satisfied
		double curL = 0;
		double curP = 0;
		for (i = 0; i < NUM_SITES; i++)
		{
			curL += f.pts[i].r;
			curP += f.pts[i].p;
		}
		double deltaL = (L - curL) / NUM_SITES;
		double deltaP = -curP / NUM_SITES;
		for (i = 0; i < NUM_SITES; i++)
		{
			f.pts[i].r += deltaL;
			f.pts[i].p += deltaP;
		}
		// check
		curL = 0;
		curP = 0;
		for (i = 0; i < NUM_SITES; i++)
		{
			curL += f.pts[i].r;
			curP += f.pts[i].p;
		}
		printf("curL: %g, should be equal to %g\n", curL, L);
		printf("curP: %g, should be zero\n", curP);
	}

	// save elongations and momenta to disk
	{
		double *r = malloc(NUM_SITES * sizeof(double));
		double *p = malloc(NUM_SITES * sizeof(double));
		for (i = 0; i < NUM_SITES; i++)
		{
			r[i] = f.pts[i].r;
			p[i] = f.pts[i].p;
		}

		WriteData("../test/collision_test_field_r0.dat", r, sizeof(double), NUM_SITES, false);
		WriteData("../test/collision_test_field_p0.dat", p, sizeof(double), NUM_SITES, false);

		free(p);
		free(r);
	}

	// event table
	event_table_t tab;
	assert(NUM_SLOTS == 16);
	EventTable_Allocate(0.1, &tab);		// dt, table

	// initial kinetic and potential energy
	double en0 = Field_Energy(&f);
	printf("initial energy: %g\n", en0);

	// collect collision events
	CollectEvents(&f, &tab);

	// process events
	SweepTable(&f, &tab);

	// final kinetic and potential energy
	double en1 = Field_Energy(&f);
	printf("final energy: %g\n", en1);
	printf("energy difference: %g, should be zero\n", en1 - en0);

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
		evtype[i] = (int)node->ev.type;
		node = node->next;
		i++;
	}
	assert(i == event_log.num);

	// save event data to disk
	WriteData("../test/collision_test_event_time.dat", evtime, sizeof(double), event_log.num, false);
	WriteData("../test/collision_test_event_iprt.dat", eviprt, sizeof(int),    event_log.num, false);
	WriteData("../test/collision_test_event_type.dat", evtype, sizeof(int),    event_log.num, false);

	// read reference data from disk
	double *evtime_ref = malloc(event_log.num * sizeof(double));
	int    *eviprt_ref = malloc(event_log.num * sizeof(int));
	int    *evtype_ref = malloc(event_log.num * sizeof(int));
	ReadData("../test/collision_test_ref_event_time.dat", evtime_ref, sizeof(double), event_log.num);
	ReadData("../test/collision_test_ref_event_iprt.dat", eviprt_ref, sizeof(int),    event_log.num);
	ReadData("../test/collision_test_ref_event_type.dat", evtype_ref, sizeof(int),    event_log.num);

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
