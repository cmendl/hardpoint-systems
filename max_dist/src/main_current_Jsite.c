/// \file main_current_Jsite.c
/// \brief Main simulation file for averaging spatial correlations of site-specific currents
/// and estimating the distribution of time-integrated currents.
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
	unsigned int ic, is, im, ip, it;

	const double pressure = 0.125;		// pressure
	const double beta = 2.0;			// inverse temperature

	const double mass[2] = { 1, 3 };	// alternating masses

	const double c_sound = 1.74374053186827856;		// speed of sound for these concrete pressure, beta and mass values

	const double ell_avr = AverageElongation(pressure, beta);
	const double  en_avr = AverageEnergy(beta);

	const double matR_0r =  0.35173188086453352255175;	// the (0,r) entry of the R matrix (normal mode 0 and physical coordinate elongation 'r')
	const double matR_1r = -2.440677042738814066293;	// the (1,r) entry of the R matrix (normal mode 1 and physical coordinate elongation 'r')
	const double matR_1p =  0.7071067811865475244008;	// the (1,p) entry of the R matrix (normal mode 1 and physical coordinate momentum 'p'): 1/sqrt(2)
	const double matR_0e =  2.8138550469162681804140;	// the (0,e) entry of the R matrix (normal mode 0 and physical coordinate energy 'e')
	const double matR_1e =  0.20275573351183709160759;	// the (1,e) entry of the R matrix (normal mode 1 and physical coordinate energy 'e')

	// number of 'sweeps' of the event table, i.e., current is integrated up to t = nsweeps*period where 'period' is the event table time period
	const unsigned int nsweeps = 1024;

	// inverse bin width for time-integrated current
	const double inv_bin_width_J = 2;
	// number of bins for time-integrated current
	#define NUM_JINT_BINS 2048
	#define MASK_JINT_BIN (NUM_JINT_BINS - 1)

	printf("pressure: %g\n", pressure);
	printf("beta:     %g\n", beta);

	printf("average elongation: %g\n", ell_avr);
	printf("average energy:     %g\n",  en_avr);

	printf("alternating masses: (%g, %g)\n", mass[0], mass[1]);

	printf("N: %i\n", NUM_SITES);

	printf("maximum time in units of event table time period: %i\n", nsweeps);

	printf("bin width for time-integrated current: %g, number of bins: %i\n", 1/inv_bin_width_J, NUM_JINT_BINS);

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
	sprintf(path, "../output/current_Jsite/sim_%i", itime);
	while (makedir(path) < 0)
	{
		printf("cannot create output directory '%s', changing 'itime'...\n", path);
		#ifdef _WIN32
		Sleep(15000);
		#else
		sleep(15);
		#endif
		itime += (clock() % 16) + 15;
		sprintf(path, "../output/current_Jsite/sim_%i", itime);
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

	const double t_max = nsweeps * period;
	printf("maximum simulation time: %g\n", t_max);

	// c_sound * t
	const double ct = c_sound * t_max;
	// rounded version
	const unsigned int cti = (unsigned int)round(ct);

	// spatial correlation structure
	spatial_correlation_t scorr;
	SpatialCorrelation_Allocate(&scorr);

	current_t *J_int = fftw_malloc(NUM_SITES * sizeof(current_t));
	currnrm_t *J_inr = fftw_malloc(NUM_SITES * sizeof(currnrm_t));
	current_t *F_int = fftw_malloc(NUM_SITES * sizeof(current_t));

	// start timer
	clock_t t_start = clock();

	for (ic = 0; ic < ncycles; ic++)
	{
		printf("cycle %i\n", ic + 1);

		// sort time-integrated current values transformed to normal modes into bins
		unsigned long long J_int_counts[3][NUM_JINT_BINS] = { 0 };

		for (is = 0; is < nsamples; is++)
		{
			// total time-integrated currents in physical coordinates for each lattice site
			memset(J_int, 0, NUM_SITES * sizeof(current_t));

			// initialize elongations and momenta using single-site probability density, Eq. (2.11)
			EquilibrateCanonical(pressure, beta, &seed, &f);

			// initial energy
			double en0 = Field_Energy(&f);
			//printf("initial energy: %g\n", en0);

			// spatial integrals over field variables with length c*t
			memset(F_int, 0, NUM_SITES * sizeof(current_t));
			for (ip = 0; ip < cti; ip++)
			{
				F_int[0].r +=                                   f.pts[(0 - ip) & SITE_MASK].r - ell_avr;
				F_int[0].p +=                  f.mass[ip & 1] * f.pts[(0 - ip) & SITE_MASK].v;
				F_int[0].e += LattPoint_Energy(f.mass[ip & 1], &f.pts[(0 - ip) & SITE_MASK])  -  en_avr;
			}
			for (ip = 1; ip < NUM_SITES; ip++)
			{
				F_int[ip].r = F_int[ip - 1].r +                                   f.pts[ip].r -                                           f.pts[(ip - cti) & SITE_MASK].r;
				F_int[ip].p = F_int[ip - 1].p +                  f.mass[ip & 1] * f.pts[ip].v -                  f.mass[(ip - cti) & 1] * f.pts[(ip - cti) & SITE_MASK].v;
				F_int[ip].e = F_int[ip - 1].e + LattPoint_Energy(f.mass[ip & 1], &f.pts[ip])  - LattPoint_Energy(f.mass[(ip - cti) & 1], &f.pts[(ip - cti) & SITE_MASK]);
			}

			// collect initial collision events
			CollectEvents(&f, &tab);

			// run actual simulation
			for (it = 0; it < nsweeps; it++)
			{
				// process events and record current
				current_t J[NUM_SITES];
				SweepTableJSite(&f, &tab, J);

				// integrate currents over time; subtract theoretical average value
				for (ip = 0; ip < NUM_SITES; ip++)
				{
					J_int[ip].r += J[ip].r;
					J_int[ip].p += J[ip].p - period*pressure;
					J_int[ip].e += J[ip].e;
				}

				// subtract 'period' from all timing events
				GlobalTimeshift(period, &f, &tab);
				// correspondingly add 'period' to time offset of event log
				#ifdef _DEBUG
				event_log.time_offset += period;
				#endif

			}	// time point loop

			// complete velocity integration for elongation current; 'period' has already been subtracted
			for (ip = 0; ip < NUM_SITES; ip++)
			{
				lattpoint_t *point = &f.pts[ip];
				J_int[ip].r -= (0 - point->t) * point->v;
			}

			// not required here
			// Field_Synchronize(0, &f);

			// transform integrated current values to normal modes
			for (ip = 0; ip < NUM_SITES; ip++)
			{
				// transform to normal modes -1, 0, 1
				J_inr[ip].m[0] = matR_1r * J_int[ip].r - matR_1p * J_int[ip].p + matR_1e * J_int[ip].e;		// mode -1
				J_inr[ip].m[1] = matR_0r * J_int[ip].r                         + matR_0e * J_int[ip].e;		// mode  0
				J_inr[ip].m[2] = matR_1r * J_int[ip].r + matR_1p * J_int[ip].p + matR_1e * J_int[ip].e;		// mode  1
			}

			// collect and accumulate spatial time-integrated current auto-correlations
			SpatialCorrelation_Collect((pointdata_t *)J_inr, (pointdata_t *)J_inr, &scorr);

			// subtract spatial integrals over field variables at t = 0 from sound modes and sort into bins
			for (ip = 0; ip < NUM_SITES; ip++)
			{
				// correspondingly transform integrals over field variables
				unsigned int iq = (ip + cti) & SITE_MASK;	// integrate "backwards"
				double Fm1 = matR_1r * F_int[iq].r - matR_1p * F_int[iq].p + matR_1e * F_int[iq].e;
				double F_1 = matR_1r * F_int[ip].r + matR_1p * F_int[ip].p + matR_1e * F_int[ip].e;

				// subtract spatial integrals over field variables at t = 0 from sound modes
				J_inr[ip].m[0] += Fm1;		// mode -1
				J_inr[ip].m[2] -= F_1;		// mode  1

				for (im = 0; im < 3; im++)	// for each mode...
				{
					// bit masking modulo operation works for negative integers also
					int ib = (int)floor(J_inr[ip].m[im] * inv_bin_width_J);
					if (-NUM_JINT_BINS/2 <= ib && ib < NUM_JINT_BINS/2)
					{
						J_int_counts[im][ib & MASK_JINT_BIN]++;
					}
				}
			}

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
			printf("average time difference between events: %g\n", nsweeps * period / event_log.num);
			// reset event log
			EventLog_Delete(&event_log);
			#endif

		}	// sample loop

		// divide by 'nsamples'
		SpatialCorrelation_ScalarMultiply(1.0 / nsamples, &scorr);

		// save to disk
		printf("saving histograms and spatial correlations of time-integrated currents to disk...\n");
		for (im = 0; im < 3; im++)	// for each mode...
		{
			const char *mode_names[3] = { "N", "0", "1" };
			sprintf(path, "../output/current_Jsite/sim_%i/sim_%i_ns%i_J_int_mode%s_hist_%04d.dat", itime, itime, nsamples, mode_names[im], ic);
			WriteData(path, J_int_counts[im], sizeof(unsigned long long), NUM_JINT_BINS, false);
		}
		// scorr.avr[1] same as scorr.avr[0]
		sprintf(path, "../output/current_Jsite/sim_%i/sim_%i_ns%i_J_int_corr_avr_%04d.dat", itime, itime, nsamples, ic); WriteData(path, &scorr.avr[0], sizeof(pointdata_t), 1, false);
		sprintf(path, "../output/current_Jsite/sim_%i/sim_%i_ns%i_J_int_corr_sqr_%04d.dat", itime, itime, nsamples, ic); WriteData(path, scorr.sqr, sizeof(cov_mat_t), NUM_SITES, false);

		// reset spatial correlation structure
		SpatialCorrelation_Reset(&scorr);

	}	// cycle loop

	clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	printf("finished simulation, CPU time: %g, average CPU time per run: %g\n", cpu_time, cpu_time / (ncycles * nsamples));

	// clean up
	fftw_free(J_inr);
	fftw_free(F_int);
	fftw_free(J_int);
	SpatialCorrelation_Delete(&scorr);
	EventTable_Delete(&tab);
	fftw_cleanup();

	return 0;
}
