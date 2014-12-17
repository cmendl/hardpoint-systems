/// \file main_fields_micro.c
/// \brief Main simulation file for averaging spatial correlations of the field variables,
/// initializing field variables according to microcanonical ensemble statistics.
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

// for sleep function and creating directories
#ifdef _WIN32
#include <windows.h>
#include <direct.h>

static int makedir(const char *path)
{
	return _mkdir(path);
}

#else

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static int makedir(const char *path)
{
	return mkdir(path, 0755);
}

#endif

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif



int main(int argc, const char* argv[])
{
	unsigned int ic, is, ip, im, it;

	const double length = 1;			// elongation (or length, stretch)
	const double energy = 1;			// average energy

	const double mass[2] = { 1, 3 };	// alternating masses

	#define NUM_SPATIAL_CORR 3			// number of time points at which spatial correlations are collected

	// cumulative number of 'sweeps' of the event table
	const unsigned int nsweeps[NUM_SPATIAL_CORR] = { 256, 512, 1024 };

	printf("length: %g\n", length);
	printf("energy: %g\n", energy);

	printf("alternating masses: (%g, %g)\n", mass[0], mass[1]);

	printf("N: %i\n", NUM_SITES);

	printf("number of sweeps of the event table: %i\n", nsweeps[NUM_SPATIAL_CORR-1]);

	// check number of input arguments
	if (argc != 3) {
		fprintf(stderr, "syntax: %s <number of cycles> <number of samples per cycle>\n", argv[0]);
		return -1;
	}

	// number of cycles
	unsigned int ncycles = atoi(argv[1]);
	printf("number of cycles: %i\n", ncycles);

	// number of samples per cycle
	unsigned int nsamples = atoi(argv[2]);
	printf("number of samples per cycle: %i\n", nsamples);

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// UNIX time
	unsigned int itime = (unsigned int)time(NULL);
	printf("itime: %i\n", itime);

	// trying to create output directory corresponding to 'itime'
	char path[1024];
	sprintf(path, "../output/fields_micro/sim_%i", itime);
	while (makedir(path) < 0)
	{
		printf("cannot create output directory '%s', changing 'itime'...\n", path);
		#ifdef _WIN32
		Sleep(15000);
		#else
		sleep(15);
		#endif
		itime += (clock() % 16) + 15;
		sprintf(path, "../output/fields_micro/sim_%i", itime);
	}
	printf("created output directory '%s'...\n", path);

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// lattice field
	field_t f;
	Field_Create(mass, &f);

	// event table
	event_table_t tab;
	EventTable_Allocate(1.0/NUM_SLOTS, &tab);	// dt, table
	printf("event table number of slots: %i, slot interval dt: %g, events per slot: %i\n", NUM_SLOTS, tab.dt, EVENTS_PER_SLOT);
	// time period
	const double period = NUM_SLOTS * tab.dt;
	printf("event table time period: %g\n", period);

	printf("maximum simulation time: %g\n", nsweeps[NUM_SPATIAL_CORR-1] * period);

	// create subdirectories for saving files
	for (im = 0; im < NUM_SPATIAL_CORR; im++)
	{
		sprintf(path, "../output/fields_micro/sim_%i/t%i", itime, (int)(nsweeps[im]*period));
		if (makedir(path) < 0) {
			fprintf(stderr, "'mkdir(%s)' failed.\n", path);
			return -1;
		}
	}

	// spatial correlation structure
	spatial_correlation_t scorr[NUM_SPATIAL_CORR];
	for (im = 0; im < NUM_SPATIAL_CORR; im++) {
		SpatialCorrelation_Allocate(&scorr[im]);
	}
	// lattice point data
	pointdata_t data_ref[NUM_SITES];
	pointdata_t data_cur[NUM_SITES];

	// start timer
	clock_t t_start = clock();

	for (ic = 0; ic < ncycles; ic++)
	{
		printf("cycle %i\n", ic + 1);

		for (is = 0; is < nsamples; is++)
		{
			// initialize elongations and momenta according to microcanonical ensemble
			EquilibrateMicrocanonical(length, energy, &seed, &f);

			// initial energy
			double en0 = Field_Energy(&f);
			//printf("initial energy: %g\n", en0);

			// initial field: store elongation and momentum, and calculate energies of field variables
			for (ip = 0; ip < NUM_SITES; ip++)
			{
				const lattpoint_t *point = &f.pts[ip];
				assert(point->t == 0);
				data_ref[ip].r = point->r;								// elongation
				data_ref[ip].p = mass[ip & 1] * point->v;				// momentum
				data_ref[ip].e = LattPoint_Energy(mass[ip & 1], point);	// energy
			}

			// collect initial collision events
			CollectEvents(&f, &tab);

			// run actual simulation
			for (im = 0; im < NUM_SPATIAL_CORR; im++)
			{
				for (it = (im == 0 ? 0 : nsweeps[im-1]); it < nsweeps[im]; it++)	// repeat several sweeps
				{
					// process events
					SweepTable(&f, &tab);

					// synchronize up to period of event table
					//Field_Synchronize(period, &f);

					// subtract 'period' from all timing events
					GlobalTimeshift(period, &f, &tab);
					// correspondingly add 'period' to time offset of event log
					#ifdef _DEBUG
					event_log.time_offset += period;
					#endif
				}

				// 'period' has already been subtracted
				Field_Synchronize(0, &f);

				// current field: store elongation and momentum, and calculate energies of field variables
				for (ip = 0; ip < NUM_SITES; ip++)
				{
					const lattpoint_t *point = &f.pts[ip];
					assert(point->t == 0);
					data_cur[ip].r = point->r;								// elongation
					data_cur[ip].p = mass[ip & 1] * point->v;				// momentum
					data_cur[ip].e = LattPoint_Energy(mass[ip & 1], point);	// energy
				}

				// collect and accumulate spatial correlations (for time offset 'nsweeps[im] * period')
				SpatialCorrelation_Collect(data_ref, data_cur, &scorr[im]);

			}	// time point loop

			// final energy
			double en1 = Field_Energy(&f);
			//printf("final energy: %g\n", en1);
			//printf("energy difference: %g (should be zero)\n", en1 - en0);
			if (fabs(en1 - en0) > 1e-10) {
				fprintf(stderr, "warning: energy seems not to be conserved, difference: %g\n", fabs(en1 - en0));
			}

			// reset event table
			EventTable_Reset(&tab);

			#ifdef _DEBUG
			printf("total number of events: %i\n", event_log.num);
			printf("average time difference between events: %g\n", nsweeps[NUM_SPATIAL_CORR-1] * period / event_log.num);
			// reset event log
			EventLog_Delete(&event_log);
			#endif

		}	// sample loop

		for (im = 0; im < NUM_SPATIAL_CORR; im++)
		{
			// divide by 'nsamples'
			SpatialCorrelation_ScalarMultiply(1.0 / nsamples, &scorr[im]);

			// save to disk
			const int timepoint = (int)(nsweeps[im]*period);
			printf("saving averages required for correlators at time %i to disk...\n", timepoint);
			sprintf(path, "../output/fields_micro/sim_%i/t%i/sim_%i_t%i_ns%i_avr0_%04d.dat", itime, timepoint, itime, timepoint, nsamples, ic); WriteData(path, &scorr[im].avr[0], sizeof(pointdata_t), 1, false);
			sprintf(path, "../output/fields_micro/sim_%i/t%i/sim_%i_t%i_ns%i_avr1_%04d.dat", itime, timepoint, itime, timepoint, nsamples, ic); WriteData(path, &scorr[im].avr[1], sizeof(pointdata_t), 1, false);
			sprintf(path, "../output/fields_micro/sim_%i/t%i/sim_%i_t%i_ns%i_sqrs_%04d.dat", itime, timepoint, itime, timepoint, nsamples, ic); WriteData(path, scorr[im].sqr, sizeof(cov_mat_t), NUM_SITES, false);

			// reset spatial correlation structure
			SpatialCorrelation_Reset(&scorr[im]);
		}

	}	// cycle loop

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("finished simulation, CPU time: %g, average CPU time per run: %g\n", cpu_time, cpu_time / (ncycles * nsamples));

	// clean up
	for (im = 0; im < NUM_SPATIAL_CORR; im++) {
		SpatialCorrelation_Delete(&scorr[im]);
	}
	EventTable_Delete(&tab);
	fftw_cleanup();

	return 0;
}
