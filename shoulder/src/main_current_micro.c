/// \file main_current_micro.c
/// \brief Main simulation file for averaging current time-correlations,
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
#include "time_correlation.h"
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
	unsigned int ic, is, im, it;

	const double length = 1.2456883086063404313235;		// elongation (or length, stretch)
	const double energy = 0.4889613325006079355656;		// average energy

	const double pressure = 1.2;						// corresponding pressure

	// maximum current correlation time, in units of event table time period
	const unsigned int n_max_corr = 1024;

	// number of event table time periods to additionally average the current over
	#define NUM_AVR_CURRENT 8

	printf("length: %g\n", length);
	printf("energy: %g\n", energy);

	printf("N: %i\n", NUM_SITES);

	printf("maximum correlation time in units of event table time period: %i\n", n_max_corr);
	printf("current is averaged over 1 and %i event table time periods\n", NUM_AVR_CURRENT);

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
	sprintf(path, "../output/current_micro/sim_%i", itime);
	while (makedir(path) < 0)
	{
		printf("cannot create output directory '%s', changing 'itime'...\n", path);
		#ifdef _WIN32
		Sleep(15000);
		#else
		sleep(15);
		#endif
		itime += (clock() % 16) + 15;
		sprintf(path, "../output/current_micro/sim_%i", itime);
	}
	printf("created output directory '%s'...\n", path);

	// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
	randseed_t seed;
	Random_SeedInit(1865811235122147685LL * (uint64_t)itime, &seed);

	// lattice field
	field_t f;
	Field_Create(&f);

	// event table
	event_table_t tab;
	EventTable_Allocate(1.0/NUM_SLOTS, &tab);	// dt, table
	printf("event table number of slots: %i, slot interval dt: %g, events per slot: %i\n", NUM_SLOTS, tab.dt, EVENTS_PER_SLOT);
	// time period
	const double period = NUM_SLOTS * tab.dt;
	printf("event table time period: %g\n", period);

	printf("maximum correlation time: %g\n", n_max_corr * period);

	// currents for each time point
	current_t *Jt = fftw_malloc(2*n_max_corr * sizeof(current_t));
	// time correlation structure for calculating current correlations
	time_correlation_t J_tcorr;
	TimeCorrelation_Allocate(n_max_corr, &J_tcorr);

	// currents additionally averaged over 'NUM_AVR_CURRENT' event table time periods
	current_t *Jt_longdt = fftw_malloc(2*n_max_corr/NUM_AVR_CURRENT * sizeof(current_t));
	// time correlation structure with 'NUM_AVR_CURRENT * period' resolution for calculating current correlations
	time_correlation_t J_tcorr_longdt;
	TimeCorrelation_Allocate(n_max_corr/NUM_AVR_CURRENT, &J_tcorr_longdt);

	// total time-integrated current
	current_t *J_int = fftw_malloc(2*nsamples * sizeof(current_t));

	// start timer
	clock_t t_start = clock();

	for (ic = 0; ic < ncycles; ic++)
	{
		printf("cycle %i\n", ic + 1);

		// initialize with zeros
		memset(J_int, 0, 2*nsamples * sizeof(current_t));

		for (is = 0; is < nsamples; is++)
		{
			// initialize elongations and momenta according to microcanonical ensemble
			EquilibrateMicrocanonical(length, energy, &seed, &f);

			// initial energy
			double en0 = Field_Energy(&f);
			//printf("initial energy: %g\n", en0);

			// collect initial collision events
			CollectEvents(&f, &tab);

			// run actual simulation
			for (it = 0; it < 2*n_max_corr; it++)	// twice the correlation time interval to collect all time differences equally often
			{
				// process events and record current
				Jt[it] = SweepTable(&f, &tab);
				// synchronize up to period of event table to complete elongation current integration
				current_t J_add = Field_Synchronize(period, &f);
				Jt[it].r += J_add.r;
				Jt[it].p += J_add.p;
				Jt[it].e += J_add.e;

				// integrate currents over time; subtract theoretical average value
				im = 2*is + (it/n_max_corr);	// separate into first and second correlation time interval
				J_int[im].r += Jt[it].r;
				J_int[im].p += Jt[it].p - pressure;
				J_int[im].e += Jt[it].e;

				// subtract 'period' from all timing events
				GlobalTimeshift(period, &f, &tab);
				// correspondingly add 'period' to time offset of event log
				#ifdef _DEBUG
				event_log.time_offset += period;
				#endif

			}	// time point loop

			// current correlations
			TimeCorrelation_Collect((pointdata_t *)Jt, &J_tcorr);

			// averaged current correlations
			for (im = 0; im < 2*J_tcorr_longdt.num_t; im++)
			{
				Jt_longdt[im].r = 0;
				Jt_longdt[im].p = 0;
				Jt_longdt[im].e = 0;

				for (it = 0; it < NUM_AVR_CURRENT; it++)
				{
					Jt_longdt[im].r += Jt[im*NUM_AVR_CURRENT + it].r;
					Jt_longdt[im].p += Jt[im*NUM_AVR_CURRENT + it].p;
					Jt_longdt[im].e += Jt[im*NUM_AVR_CURRENT + it].e;
				}

				// normalization
				Jt_longdt[im].r *= 1.0 / NUM_AVR_CURRENT;
				Jt_longdt[im].p *= 1.0 / NUM_AVR_CURRENT;
				Jt_longdt[im].e *= 1.0 / NUM_AVR_CURRENT;
			}
			TimeCorrelation_Collect((pointdata_t *)Jt_longdt, &J_tcorr_longdt);

			// final energy
			double en1 = Field_Energy(&f);
			//printf("final energy: %g\n", en1);
			//printf("energy difference: %g (should be zero)\n", en1 - en0);
			if (fabs(en1 - en0) > 1e-10) {
				fprintf(stderr, "warning: energy seems not to be conserved, difference: %g\n", fabs(en1 - en0));
			}

			// SweepTable() returns current per time
			J_int[2*is  ].r *= period;
			J_int[2*is  ].p *= period;
			J_int[2*is  ].e *= period;
			J_int[2*is+1].r *= period;
			J_int[2*is+1].p *= period;
			J_int[2*is+1].e *= period;

			// reset event table
			EventTable_Reset(&tab);

			#ifdef _DEBUG
			// reset event log
			EventLog_Delete(&event_log);
			#endif

		}	// sample loop

		// divide by 'nsamples'
		TimeCorrelation_ScalarMultiply(1.0 / nsamples, &J_tcorr);
		TimeCorrelation_ScalarMultiply(1.0 / nsamples, &J_tcorr_longdt);

		// save to disk
		printf("saving average currents and current correlations to disk...\n");
		sprintf(path, "../output/current_micro/sim_%i/sim_%i_ns%i_Javr_%04d.dat", itime, itime, nsamples, ic); WriteData(path, J_tcorr.avr, sizeof(current_t), 2*n_max_corr, false);
		sprintf(path, "../output/current_micro/sim_%i/sim_%i_ns%i_Jsqr_%04d.dat", itime, itime, nsamples, ic); WriteData(path, J_tcorr.sqr, sizeof(cov_mat_t),   n_max_corr, false);

		sprintf(path, "../output/current_micro/sim_%i/sim_%i_ns%i_dt%d_Javr_%04d.dat", itime, itime, nsamples, NUM_AVR_CURRENT, ic); WriteData(path, J_tcorr_longdt.avr, sizeof(current_t), 2*J_tcorr_longdt.num_t, false);
		sprintf(path, "../output/current_micro/sim_%i/sim_%i_ns%i_dt%d_Jsqr_%04d.dat", itime, itime, nsamples, NUM_AVR_CURRENT, ic); WriteData(path, J_tcorr_longdt.sqr, sizeof(cov_mat_t),   J_tcorr_longdt.num_t, false);

		printf("saving integrated current to disk...\n");
		sprintf(path, "../output/current_micro/sim_%i/sim_%i_ns%i_Jint_%04d.dat", itime, itime, nsamples, ic); WriteData(path, J_int, sizeof(current_t), 2*nsamples, false);

		// reset time correlation structures
		TimeCorrelation_Reset(&J_tcorr);
		TimeCorrelation_Reset(&J_tcorr_longdt);

	}	// cycle loop

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("finished simulation, CPU time: %g, average CPU time per run: %g\n", cpu_time, cpu_time / (ncycles * nsamples));

	// clean up
	fftw_free(J_int);
	TimeCorrelation_Delete(&J_tcorr_longdt);
	TimeCorrelation_Delete(&J_tcorr);
	fftw_free(Jt_longdt);
	fftw_free(Jt);
	EventTable_Delete(&tab);
	fftw_cleanup();

	return 0;
}
