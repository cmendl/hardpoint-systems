/// \file field.c
/// \brief Lattice point and lattice field data structures.
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

#include "field.h"
#include "potential.h"
#include "util.h"
#include <stdlib.h>
#include <memory.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Create field (initialize masses, set elongations and momenta to zero)
///
void Field_Create(const double mass[2], field_t *f)
{
	unsigned int i;

	// initialize with zeros
	memset(f->pts, 0, sizeof(f->pts));

	// initially, event indices are invalid
	for (i = 0; i < NUM_SITES; i++) {
		f->pts[i].ievent = -1;
	}

	f->mass[0] = mass[0];
	f->mass[1] = mass[1];

	// average mass
	double m_mean = 0.5*(mass[0] + mass[1]);

	f->relm[0] = mass[0] / m_mean;
	f->relm[1] = mass[1] / m_mean;

	// product of masses
	f->rm_prod = mass[0]*mass[1] / m_mean;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Advance particle elongations to time 'time' and record elongation current, assuming that no collisions occur in the meantime
///
currsum_even_odd_t Field_Synchronize(const double time, field_t *f)
{
	unsigned int i;
	currsum_even_odd_t J = { 0 };

	for (i = 0; i < NUM_SITES; i++)
	{
		// reference to lattice point
		lattpoint_t *point = &f->pts[i];

		double dt = time - point->t;
		//assert(dt >= 0);	// 'time' must be in the future (but particle position is advanced by a small time step after collision, such that point->t can be slightly larger)

		// elongation current is negative time-integrated velocity up to 'ev->time'
		J.evenodd[i & 1].r -= dt * point->v;

		// update elongation and time
		point->r += dt * (f->pts[(i+1) & SITE_MASK].v - point->v);
		point->t = time;
	}
	// normalize total elongation current
	J.evn.r /= NUM_SITES;
	J.odd.r /= NUM_SITES;

	return J;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the total elongation of the field
///
double Field_Elongation(const field_t *f)
{
	unsigned int i;

	double r = 0;
	for (i = 0; i < NUM_SITES; i++)
	{
		r += f->pts[i].r;
	}

	return r;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the total momentum of the field
///
double Field_Momentum(const field_t *f)
{
	unsigned int i;

	double p = 0;
	for (i = 0; i < NUM_SITES; i += 2)
	{
		p += f->mass[0] * f->pts[i  ].v;
		p += f->mass[1] * f->pts[i+1].v;
	}

	return p;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the total energy of a lattice point
///
double LattPoint_Energy(const double mass, const lattpoint_t *point)
{
	// kinetic energy
	double en = 0.5*mass*square(point->v);

	// potential energy is always zero (hard-core potential)
	assert(point->r >= V_WIDTH_CORE);

	return en;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the total energy of the field
///
double Field_Energy(const field_t *f)
{
	unsigned int i;

	double en = 0;
	for (i = 0; i < NUM_SITES; i += 2)
	{
		en += LattPoint_Energy(f->mass[0], &f->pts[i  ]);
		en += LattPoint_Energy(f->mass[1], &f->pts[i+1]);
	}

	return en;
}
