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
/// \brief Create field (set elongations and momenta to zero)
///
void Field_Create(field_t *f)
{
	unsigned int i;

	// initialize with zeros
	memset(f->pts, 0, sizeof(f->pts));

	// initially, event indices are invalid
	for (i = 0; i < NUM_SITES; i++) {
		f->pts[i].ievent = -1;
	}
}


//_______________________________________________________________________________________________________________________
///
/// \brief Advance particle elongations to time 'time' and record elongation current, assuming that no collisions occur in the meantime
///
current_t Field_Synchronize(const double time, field_t *f)
{
	unsigned int i;
	current_t J = { 0 };

	for (i = 0; i < NUM_SITES; i++)
	{
		// reference to lattice point
		lattpoint_t *point = &f->pts[i];

		double dt = time - point->t;
		//assert(dt >= 0);	// 'time' must be in the future (but particle position is advanced by a small time step after collision, such that point->t can be slightly larger)

		// elongation current is negative time-integrated momentum up to 'ev->time'
		J.r -= dt * point->p;

		// update elongation and time
		point->r += dt * (f->pts[(i+1) & SITE_MASK].p - point->p);
		point->t = time;
	}
	// normalize total elongation current
	J.r /= NUM_SITES;

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
	for (i = 0; i < NUM_SITES; i++)
	{
		p += f->pts[i].p;
	}

	return p;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Calculate the total energy of a lattice point
///
double LattPoint_Energy(const lattpoint_t *point)
{
	double en;

	// kinetic energy
	en = 0.5*square(point->p);

	// shoulder potential energy
	if (point->r <= V_WIDTH_PLAT)
	{
		assert(point->r >= V_WIDTH_CORE);
		en += V_POT_HEIGHT;
	}

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

	for (i = 0; i < NUM_SITES; i++)
	{
		en += LattPoint_Energy(&f->pts[i]);
	}

	return en;
}
