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

	return V_WIDTH_CORE - params->inv_q * log(1 - t);
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average elongation (or length) for given pressure and inverse temperature beta
///
double AverageElongation(const double pressure, const double beta)
{
	const double q = pressure * beta;

	return V_WIDTH_CORE + 1 / q;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average energy for given pressure and inverse temperature beta
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


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average pressure for given elongation and energy
///
double AveragePressure(const double elongation, const double energy)
{
	return 2*energy/elongation;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the average inverse temperature beta for given energy
///
double AverageBeta(const double energy)
{
	return 0.5/energy;
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

	// 1 / m_i
	const double inv_mass[2] = { 1/f->mass[0], 1/f->mass[1] };

	unsigned int i, j, k, l, n;

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
			double p_avr = sqrt((en_loc[i] + en_loc[i+1]) * f->rm_prod);
			double v0 =  p_avr * inv_mass[0];
			double v1 = -p_avr * inv_mass[1];	// total momentum is zero

			f->pts[i  ].v = v0;
			f->pts[i+1].v = v1;
		}
	}

	// equilibrate elongations (independent of momenta)
	for (n = 0; n < NUM_PASSES; n++)
	{
		for (i = 0; i < NUM_SITES; i++)
		{
			// choose partner randomly, must be different from i
			while((j = Random_GetUint(seed) & SITE_MASK) == i)
				;

			// local sum is conserved
			double r_sum = f->pts[i].r + f->pts[j].r - 2*V_WIDTH_CORE;
			double dri = r_sum * Random_GetUniform(seed);
			f->pts[i].r = V_WIDTH_CORE + dri;
			f->pts[j].r = V_WIDTH_CORE + (r_sum - dri);

			// must be non-negative
			assert(f->pts[i].r >= V_WIDTH_CORE);
			assert(f->pts[j].r >= V_WIDTH_CORE);
		}
	}
	// check field elongation
	assert(fabs(Field_Elongation(f)/NUM_SITES - elongation) < 1e-12);

	// precompute expressions depending on masses
	const double sqrt_mass[2] = { sqrt(2*f->mass[0]), sqrt(2*f->mass[1]) };
	const double sqrt_mtwo[2][2] = {
		{ sqrt(2*f->mass[0]), sqrt(6*f->mass[0]*f->mass[1] / (2*f->mass[0] + f->mass[1])) },
		{ sqrt(6*f->mass[0]*f->mass[1] / (f->mass[0] + 2*f->mass[1])), sqrt(2*f->mass[1]) }
	};

	// equilibrate momenta (independent of elongations)
	for (n = 0; n < NUM_PASSES; n++)
	{
		for (l = 0; l < NUM_SITES; l++)
		{
			// choose two partners randomly
			do
			{
				j = Random_GetUint(seed) & SITE_MASK;
				k = Random_GetUint(seed) & SITE_MASK;
			}
			while(l == j || l == k || j == k);	// must be pairwise different

			// permute indices such that i (== l) and k have same masses
			if (((l - k) & 1) == 0)	// l == k mod 2
			{
				// simply copy index
				i = l;
			}
			else if (((l - j) & 1) == 0)	// l == j mod 2
			{
				i = l;
				// flip j <-> k
				unsigned int t = k;
				k = j;
				j = t;
			}
			else	// j == k mod 2
			{
				assert(((j - k) & 1) == 0);
				// flip i <-> j
				i = j;
				j = l;
			}
			assert(((i - k) & 1) == 0);
			// must be pairwise different
			assert(i != j && i != k && j != k);

			// local masses (can be equal)
			const double lmass[3] = {
				f->mass[i & 1],
				f->mass[j & 1],
				f->mass[k & 1],
			};
			// corresponding inverse local masses
			const double inv_lmass[3] = {
				inv_mass[i & 1],
				inv_mass[j & 1],
				inv_mass[k & 1],
			};

			// local total momentum (must be conserved)
			const double p_tot =
				lmass[0] * f->pts[i].v +
				lmass[1] * f->pts[j].v +
				lmass[2] * f->pts[k].v;

			// local total energy (must be conserved)
			double e_tot = 0.5*(
				lmass[0] * square(f->pts[i].v) +
				lmass[1] * square(f->pts[j].v) +
				lmass[2] * square(f->pts[k].v));

			// 1 / (m_i + m_j + m_k)
			double inv_tot_mass = 1 / (lmass[0] + lmass[1] + lmass[2]);

			// center of the ellipse (intersection of energy ellipsoid with momentum plane)
			double ctr[3] = {
				p_tot * inv_tot_mass * lmass[0],
				p_tot * inv_tot_mass * lmass[1],
				p_tot * inv_tot_mass * lmass[2],
			};

			// intrinsic energy: difference between total energy and center of lmass kinetic energy
			double en_int = e_tot - 0.5*inv_tot_mass*square(p_tot);
			assert(en_int >= 0);	// must be non-negative
			en_int = sqrt(en_int);

			double lambda[2] = { sqrt_mass[i & 1], sqrt_mtwo[i & 1][j & 1] };

			// random angle phi
			double phi = M_2PI * Random_GetUniform(seed);
			double lambda_cphi = lambda[0] * cos(phi);
			double lambda_sphi = lambda[1] * sin(phi);

			// update velocities corresponding to new momenta
			f->pts[i].v = (ctr[0] + en_int*(lambda_cphi*u[0] + lambda_sphi*v[0])) * inv_lmass[0];
			f->pts[j].v = (ctr[1] + en_int*(lambda_cphi*u[1] + lambda_sphi*v[1])) * inv_lmass[1];
			f->pts[k].v = (ctr[2] + en_int*(lambda_cphi*u[2] + lambda_sphi*v[2])) * inv_lmass[2];

			// check momentum and energy conservation
			assert(fabs(lmass[0]*f->pts[i].v + lmass[1]*f->pts[j].v + lmass[2]*f->pts[k].v - p_tot) < 1e-12);
			assert(fabs(0.5*(
				lmass[0] * square(f->pts[i].v) +
				lmass[1] * square(f->pts[j].v) +
				lmass[2] * square(f->pts[k].v)) - e_tot) < 1e-12);
		}
	}
	// check field momentum and energy
	assert(fabs(Field_Momentum(f)/NUM_SITES) < 1e-12);
	assert(fabs(Field_Energy(f)/NUM_SITES - energy) < 1e-12);

	// initially, event indices are invalid
	for (i = 0; i < NUM_SITES; i++) {
		f->pts[i].ievent = -1;
	}
}
