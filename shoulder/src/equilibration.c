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
/// \brief Fill 'ensemble_params_t' structure given the pressure and inverse temperature beta
///
void EnsembleParams_Fill(const double pressure, const double beta, ensemble_params_t *params)
{
	const double q = pressure * beta;

	params->hbeta = V_POT_HEIGHT*beta;

	const double tmp = exp(-params->hbeta) * (exp((V_WIDTH_PLAT-V_WIDTH_CORE)*q) - 1);

	params->inv_q = 1 / q;
	params->tcrit = 1 / (1 + 1 / tmp);
	params->Zspat = params->inv_q * exp(-V_WIDTH_PLAT*q)*(1 + tmp);
	params->coef1 = exp(V_WIDTH_CORE*q + params->hbeta) * params->Zspat * q;
	params->coef2 = V_WIDTH_PLAT - params->inv_q * log(1 + tmp);

	params->sigma = 1 / sqrt(beta);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Evaluate the inverse cumulative distribution function (CDF)
/// of the single-site probability density for the elongation
///
double ElongationInverseCDF(const ensemble_params_t *params, const double t)
{
	assert(0 <= t && t <= 1);

	if (t < params->tcrit)
	{
		return V_WIDTH_CORE - params->inv_q * log(1 - params->coef1*t);
	}
	else
	{
		return params->coef2 - params->inv_q * log(1 - t);
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average elongation (or length) for given pressure and inverse temperature beta
///
double AverageElongation(const double pressure, const double beta)
{
	const double q = pressure * beta;

	return 1 / q + V_WIDTH_CORE + (V_WIDTH_PLAT-V_WIDTH_CORE) / (1 + exp((V_WIDTH_PLAT-V_WIDTH_CORE)*q) / (exp(V_POT_HEIGHT*beta) - 1));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average energy for given pressure and inverse temperature beta
///
double AverageEnergy(const double pressure, const double beta)
{
	const double q = pressure * beta;
	const double tmp = exp(-V_POT_HEIGHT*beta) * (exp((V_WIDTH_PLAT-V_WIDTH_CORE)*q) - 1);

	return  0.5/beta + V_POT_HEIGHT*(1 - 1 / (1 + tmp));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Initialize elongations and momenta randomly according to single-site probability density of canonical ensemble
///
void EquilibrateCanonical(const ensemble_params_t *params, randseed_t *seed, field_t *f)
{
	unsigned int i;

	memset(f->pts, 0, NUM_SITES * sizeof(lattpoint_t));

	for (i = 0; i < NUM_SITES; i++)
	{
		f->pts[i].r = ElongationInverseCDF(params, Random_GetUniform(seed));
		f->pts[i].p = params->sigma * Random_GetNormal(seed);
	}

	// initially, event indices are invalid
	for (i = 0; i < NUM_SITES; i++) {
		f->pts[i].ievent = -1;
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Initialize elongations and momenta randomly according to microcanonical ensemble,
///	for given average elongation and energy per lattice site; the total momentum is always zero
///
void EquilibrateMicrocanonical(const double elongation, const double energy, randseed_t *seed, field_t *f)
{
	// number of passes over all field variables
	#define NUM_PASSES 12

	static const double u[3] = { -0.7071067811865475, 0,                   0.7071067811865475 };	// ( -1, 0,  1 ) / sqrt(2)
	static const double v[3] = { -0.4082482904638630, 0.8164965809277260, -0.4082482904638630 };	// ( -1, 2, -1 ) / sqrt(6)

	unsigned int i, j, k, n;

	// initialize all field variables with zeros
	memset(f->pts, 0, NUM_SITES * sizeof(lattpoint_t));

	// set elongations and momenta according to microcanonical constraints
	{
		for (i = 0; i < NUM_SITES; i++)
		{
			f->pts[i].r = elongation;
		}

		// perturb local energy a bit
		double en_loc[NUM_SITES];
		double en_sum = 0;
		for (i = 0; i < NUM_SITES; i++)
		{
			en_loc[i] = energy*(0.5 + Random_GetUniform(seed));
			assert(en_loc[i] >= 0);
			en_sum += en_loc[i];
		}
		double en_diff = energy - en_sum/NUM_SITES;
		// adjust energies such that average energy equals given energy
		for (i = 0; i < NUM_SITES; i++)
		{
			en_loc[i] += en_diff;
			assert(en_loc[i] >= 0);
		}

		for (i = 0; i < NUM_SITES; i += 2)
		{
			const double p_loc = sqrt(en_loc[i] + en_loc[i+1]);

			f->pts[i  ].p =  p_loc;
			f->pts[i+1].p = -p_loc;	// total momentum is zero
		}
	}

	// perform several equilibration passes
	for (n = 0; n < NUM_PASSES; n++)
	{
		for (i = 0; i < NUM_SITES; i++)
		{
			// choose two partners randomly
			do
			{
				j = Random_GetUint(seed) & SITE_MASK;
				k = Random_GetUint(seed) & SITE_MASK;
			}
			while (i == j || i == k || j == k);	// must be pairwise different

			// local sum of relative elongations (must be conserved)
			const double r_tot = (f->pts[i].r + f->pts[j].r + f->pts[k].r) - 3*V_WIDTH_CORE;
			assert(r_tot >= 0);

			// local average momentum (must be conserved)
			const double p_avr = (f->pts[i].p + f->pts[j].p + f->pts[k].p) / 3.0;

			// local total energy (must be conserved)
			const double e_tot = 0.5*(square(f->pts[i].p) + square(f->pts[j].p) + square(f->pts[k].p)) +
				(f->pts[i].r < V_WIDTH_PLAT ? V_POT_HEIGHT : 0) +
				(f->pts[j].r < V_WIDTH_PLAT ? V_POT_HEIGHT : 0) +
				(f->pts[k].r < V_WIDTH_PLAT ? V_POT_HEIGHT : 0);

			for (; ;)	// repeat until intrinsic energy is non-negative
			{
				// equilibrate elongations
				// uniform distribution on simplex r_i + r_j + r_k = const
				double a = Random_GetUniform(seed);
				double b = Random_GetUniform(seed);
				if (1 - a - b >= 0)
				{
					f->pts[i].r = V_WIDTH_CORE + r_tot * (1 - a - b);
					f->pts[j].r = V_WIDTH_CORE + r_tot * a;
					f->pts[k].r = V_WIDTH_CORE + r_tot * b;
				}
				else
				{
					f->pts[i].r = V_WIDTH_CORE + r_tot * (a + b - 1);
					f->pts[j].r = V_WIDTH_CORE + r_tot * (1 - a);
					f->pts[k].r = V_WIDTH_CORE + r_tot * (1 - b);
				}

				// remaining kinetic energy
				const double e_kin = e_tot -
					((f->pts[i].r < V_WIDTH_PLAT ? V_POT_HEIGHT : 0) +
					 (f->pts[j].r < V_WIDTH_PLAT ? V_POT_HEIGHT : 0) +
					 (f->pts[k].r < V_WIDTH_PLAT ? V_POT_HEIGHT : 0));

				// intrinsic energy: difference between kinetic energy and center of mass kinetic energy
				const double en_int = e_kin - 3 * 0.5 * square(p_avr);
				if (en_int < 0) {
					continue;
				}

				// center of the momentum circle (intersection of kinetic energy sphere with momentum plane) is at p_avr * (1, 1, 1)

				// random angle phi
				double phi = M_2PI * Random_GetUniform(seed);
				double lambda = sqrt(2*en_int);
				double lambda_cphi = lambda * cos(phi);
				double lambda_sphi = lambda * sin(phi);

				// update momenta
				f->pts[i].p = p_avr + lambda_cphi*u[0] + lambda_sphi*v[0];
				f->pts[j].p = p_avr + lambda_cphi*u[1] + lambda_sphi*v[1];
				f->pts[k].p = p_avr + lambda_cphi*u[2] + lambda_sphi*v[2];

				// check momentum and kinetic energy conservation
				assert(fabs(f->pts[i].p + f->pts[j].p + f->pts[k].p - 3*p_avr) < 1e-12);
				assert(fabs(0.5*(square(f->pts[i].p) + square(f->pts[j].p) + square(f->pts[k].p)) - e_kin) < 1e-12);

				break;
			}
		}
	}

	// check prescribed field elongation, momentum and energy
	assert(fabs(Field_Elongation(f)/NUM_SITES - elongation) < 1e-12);
	assert(fabs(Field_Momentum(f)/NUM_SITES) < 1e-12);
	assert(fabs(Field_Energy(f)/NUM_SITES - energy) < 1e-12);

	// initially, event indices are invalid
	for (i = 0; i < NUM_SITES; i++) {
		f->pts[i].ievent = -1;
	}
}
