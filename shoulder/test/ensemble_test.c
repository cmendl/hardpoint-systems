/// \file ensemble_test.c
/// \brief Test canonical ensemble cumulative distribution function (CDF).
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
#include "util.h"
#include <stdlib.h>
#include <stdio.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


double ElongationInverseCDF(const ensemble_params_t *params, const double t);



int main()
{
	const double beta = 2.0;		// inverse temperature
	const double pressure = 1.2;	// pressure

	// reference average length and energy
	const double length_ref = 1.2456883086063404;
	const double energy_ref = 0.4889613325006078;

	// number of evaluation points
	const unsigned int num = 25;

	unsigned int i;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// average length and energy
	double length = AverageElongation(pressure, beta);
	double energy = AverageEnergy(pressure, beta);
	printf("Average length: %g, error: %g\n", length, fabs(length - length_ref));
	printf("Average energy: %g, error: %g\n", energy, fabs(energy - energy_ref));

	ensemble_params_t params;
	EnsembleParams_Fill(pressure, beta, &params);

	printf("Evaluating inverse CDF function . . .\n");

	// evaluate inverse CDF function at grid points
	double *invCDF = malloc(num * sizeof(double));
	for (i = 0; i < num; i++)
	{
		invCDF[i] = ElongationInverseCDF(&params, (double)i/num);
	}

	// save to disk
	WriteData("../test/ensemble_test.dat", invCDF, sizeof(double), num, false);

	// read reference data from disk
	double *invCDF_ref = malloc(num * sizeof(double));
	int hr = ReadData("../test/ensemble_test_ref.dat", invCDF_ref, sizeof(double), num);
	if (hr < 0) {
		return -1;
	}

	// calculate average error
	double err = 0;
	for (i = 0; i < num; i++)
	{
		err += fabs(invCDF[i] - invCDF_ref[i]);
	}
	err /= num;
	printf("Average error: %g\n", err);

	// clean up
	free(invCDF_ref);
	free(invCDF);

	return 0;
}
