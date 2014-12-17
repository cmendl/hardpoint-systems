/// \file simulation_test.c
/// \brief Test simulation.
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
#include "spatial_correlation.h"
#include "util.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


// Calculate covariance matrix
static void CalculateCovariance(const spatial_correlation_t *scorr, cov_mat_t *cov)
{
	int i, j;

	// outer product of lattice point data
	cov_mat_t prod;
	prod.entries.rr = scorr->avr[0].r * scorr->avr[1].r;
	prod.entries.pr = scorr->avr[0].p * scorr->avr[1].r;
	prod.entries.er = scorr->avr[0].e * scorr->avr[1].r;

	prod.entries.rp = scorr->avr[0].r * scorr->avr[1].p;
	prod.entries.pp = scorr->avr[0].p * scorr->avr[1].p;
	prod.entries.ep = scorr->avr[0].e * scorr->avr[1].p;

	prod.entries.re = scorr->avr[0].r * scorr->avr[1].e;
	prod.entries.pe = scorr->avr[0].p * scorr->avr[1].e;
	prod.entries.ee = scorr->avr[0].e * scorr->avr[1].e;

	for (i = 0; i < NUM_SITES; i++) {
		for (j = 0; j < 9; j++) {
			cov[i].mat[j] = scorr->sqr[i].mat[j] - prod.mat[j];
		}
	}
}


int main()
{
	unsigned int i, j;

	const double pressure = 1.2;			// pressure
	const double beta = 2.0;				// inverse temperature

	const unsigned int nsweeps = /*500*/250;	// number of 'sweeps' of the event table

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	printf("pressure: %g\n", pressure);
	printf("beta:     %g\n", beta);

	printf("N: %d\n", NUM_SITES);
	assert(NUM_SITES == 256);

	// lattice field
	field_t f;
	Field_Create(&f);
	{
		// artificial UNIX time
		unsigned int itime = 1420000241;

		// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
		randseed_t seed;
		Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

		ensemble_params_t params;
		EnsembleParams_Fill(pressure, beta, &params);

		// initialize elongations and momenta using single-site probability density, Eq. (2.11)
		EquilibrateCanonical(&params, &seed, &f);

		// current mean values
		double r_mean = 0;		// mean of r's (elongations), remains constant in time
		double p_mean = 0;		// mean momentum
		double e_mean = 0;		// mean energy
		for (i = 0; i < NUM_SITES; i++)
		{
			const lattpoint_t *point = &f.pts[i];
			r_mean += point->r;
			p_mean += point->p;
			e_mean += LattPoint_Energy(point);
		}
		r_mean /= NUM_SITES;
		p_mean /= NUM_SITES;
		e_mean /= NUM_SITES;
		//e_mean -= 0.5 * square(p_mean);	// internal energy per particle

		// compare with stochastic expectation values
		double r_avr = AverageElongation(pressure, beta);
		double p_avr = 0;
		double e_avr = AverageEnergy(pressure, beta);

		printf("elongation mean: %g, stochastic average: %g\n", r_mean, r_avr);
		printf("momentum   mean: %g, stochastic average: %g\n", p_mean, p_avr);
		printf("energy     mean: %g, stochastic average: %g\n", e_mean, e_avr);
	}

	// event table
	event_table_t tab;
	assert(NUM_SLOTS == 128/*64*/);
	EventTable_Allocate(0.01, &tab);	// dt, table
	printf("event table number of slots: %d, slot interval dt: %g\n", NUM_SLOTS, tab.dt);
	// time period
	const double period = NUM_SLOTS * tab.dt;
	printf("event table time period: %g\n", period);

	// spatial correlation structure
	spatial_correlation_t scorr;
	SpatialCorrelation_Allocate(&scorr);
	// lattice point data
	pointdata_t data_ref[NUM_SITES];
	pointdata_t data_cur[NUM_SITES];
	// actual covariance matrix
	cov_mat_t *cov = malloc(NUM_SITES * sizeof(cov_mat_t));

	// initial kinetic and potential energy
	double en0 = Field_Energy(&f);
	printf("initial energy: %g\n", en0);

	// start timer
	clock_t t_start = clock();

	// collect initial collision events
	CollectEvents(&f, &tab);

	for (i = 0; i < nsweeps; i++)	// repeat several sweeps
	{
		// store elongation and momentum, and calculate energies of field variables
		for (j = 0; j < NUM_SITES; j++)
		{
			const lattpoint_t *point = &f.pts[j];
			assert(point->t == 0);
			data_ref[j].r = point->r;					// elongation
			data_ref[j].p = point->p;					// momentum
			data_ref[j].e = LattPoint_Energy(point);	// energy
		}

		// process events
		SweepTable(&f, &tab);

		// synchronize up to period of event table
		Field_Synchronize(period, &f);

		// store elongation and momentum, and calculate energies of field variables
		for (j = 0; j < NUM_SITES; j++)
		{
			const lattpoint_t *point = &f.pts[j];
			assert(point->t == period);
			data_cur[j].r = point->r;					// elongation
			data_cur[j].p = point->p;					// momentum
			data_cur[j].e = LattPoint_Energy(point);	// energy
		}

		// collect and accumulate spatial correlations (for time offset 'period')
		SpatialCorrelation_Collect(data_ref, data_cur, &scorr);

		// subtract 'period' from all timing events
		GlobalTimeshift(period, &f, &tab);
		// correspondingly add 'period' to time offset of event log
		#ifdef _DEBUG
		event_log.time_offset += period;
		#endif
	}

	clock_t t_end = clock();
	printf("finished simulation, simulation time: %g, CPU time: %g\n", nsweeps * period, (double)(t_end - t_start) / CLOCKS_PER_SEC);

	// divide by 'nsweeps'
	SpatialCorrelation_ScalarMultiply(1.0/nsweeps, &scorr);
	// calculate covariance matrix
	CalculateCovariance(&scorr, cov);
	// save to disk
	printf("saving covariance matrix to disk...\n");
	WriteData("../test/simulation_test_cov.dat", cov, sizeof(cov_mat_t), NUM_SITES, false);

	// final kinetic and potential energy
	double en1 = Field_Energy(&f);
	printf("final energy: %g\n", en1);
	printf("energy difference: %g (should be zero)\n", en1 - en0);

	#ifdef _DEBUG
	printf("total number of events: %i\n", event_log.num);
	// save event data to disk
	i = 0;
	event_node_t *node = event_log.first;
	double *evtime = malloc(event_log.num * sizeof(double));
	int    *eviprt = malloc(event_log.num * sizeof(int));
	int    *evtype = malloc(event_log.num * sizeof(int));
	while (node != NULL)
	{
		evtime[i] = node->ev.time;
		eviprt[i] = node->ev.iprt;
		evtype[i] = (int)node->ev.type;
		node = node->next;
		i++;
	}
	assert(i == event_log.num);
	// save event data to disk
	WriteData("../test/simulation_test_event_time.dat", evtime, sizeof(double), event_log.num, false);
	WriteData("../test/simulation_test_event_iprt.dat", eviprt, sizeof(int),    event_log.num, false);
	WriteData("../test/simulation_test_event_type.dat", evtype, sizeof(int),    event_log.num, false);
	// clean up
	free(evtype);
	free(eviprt);
	free(evtime);
	EventLog_Delete(&event_log);
	#endif

	// clean up
	free(cov);
	SpatialCorrelation_Delete(&scorr);
	EventTable_Delete(&tab);
	fftw_cleanup();

	return 0;
}
