/// \file equilibration.c
/// \brief Equilibration of the lattice field variables: initialize elongations
/// and momenta randomly according to (micro-)canonical ensemble statistics.
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
#include "potential.h"
#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


#define M_2PI		6.283185307179586476925		//!< 2*pi


//_______________________________________________________________________________________________________________________
///
/// \brief Fill 'ensemble_params_t' structure given the pressure, inverse temperature beta and alternating masses
///
void EnsembleParams_Fill(const double pressure, const double beta, const double mass[2], ensemble_params_t *params)
{
	const double q = pressure * beta;

	params->inv_q = 1 / q;
	params->coeff = 1 - exp(-V_WIDTH_WELL*q);

	params->sigma[0] = 1 / sqrt(mass[0] * beta);
	params->sigma[1] = 1 / sqrt(mass[1] * beta);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Evaluate the inverse cumulative distribution function (CDF)
/// of the single-site probability density for the elongation
///
double ElongationInverseCDF(const ensemble_params_t *params, const double t)
{
	assert(0 <= t && t <= 1);

	if (isinf(params->inv_q))
	{
		return V_WIDTH_WELL * t;
	}
	else
	{
		return -params->inv_q * log(1 - params->coeff*t);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average elongation (or length) for given pressure and inverse temperature beta
///
double AverageElongation(const double pressure, const double beta)
{
	const double q = pressure * beta;

	if (q == 0)
	{
		return 0.5 * V_WIDTH_WELL;
	}
	else
	{
		return 1 / q - V_WIDTH_WELL / (exp(V_WIDTH_WELL*q) - 1);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average energy for given inverse temperature beta
///
double AverageEnergy(const double beta)
{
	return  0.5/beta;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Initialize elongations and momenta randomly according to single-site probability density of canonical ensemble
///
void EquilibrateCanonical(const double pressure, const double beta, randseed_t *seed, field_t *f)
{
	unsigned int i;

	// N must be even
	//assert(f->N % 2 == 0);

	// canonical ensemble parameters
	ensemble_params_t params;
	EnsembleParams_Fill(pressure, beta, f->mass, &params);

	memset(f->pts, 0, NUM_SITES * sizeof(lattpoint_t));

	// statistics for elongation is independent of mass
	for (i = 0; i < NUM_SITES; i++)
	{
		f->pts[i].r = ElongationInverseCDF(&params, Random_GetUniform(seed));
	}

	for (i = 0; i < NUM_SITES; i += 2)
	{
		f->pts[i  ].v = params.sigma[0] * Random_GetNormal(seed);
		f->pts[i+1].v = params.sigma[1] * Random_GetNormal(seed);
	}

	// initially, event indices are invalid
	for (i = 0; i < NUM_SITES; i++) {
		f->pts[i].ievent = -1;
	}
}
