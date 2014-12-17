//	Lattice point and lattice field data structures.
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

#ifndef FIELD_H
#define FIELD_H

#include "current.h"


//_______________________________________________________________________________________________________________________
//

#define NUM_SITES			4096				//!< N, number of lattice sites or particles
#define SITE_MASK			(NUM_SITES - 1)		//!< bit mask for fast modulo N operation; N must be a power of 2


//_______________________________________________________________________________________________________________________
///
/// \brief Lattice point
///
typedef struct
{
	double r;		//!< compression or elongation (alternatively, positional difference to next particle)
	double v;		//!< velocity, i.e., momentum p = mass * v
	double t;		//!< time stamp of position and momentum
	int ievent;		//!< event index, or -1 if no event is associated with lattice point
}
lattpoint_t;


double LattPoint_Energy(const double mass, const lattpoint_t *point);


//_______________________________________________________________________________________________________________________
///
/// \brief Lattice field
///
typedef struct
{
	lattpoint_t pts[NUM_SITES];		//!< lattice points
	double mass[2];					//!< the two alternating masses
	double relm[2];					//!< m_i / m_mean, where m_mean = (m_0 + m_1)/2
	double rm_prod;					//!< m_0 * m_1 / m_mean
}
field_t;


void Field_Create(const double mass[2], field_t *f);


//_______________________________________________________________________________________________________________________
///
/// \brief Cumulative current separately for even and odd lattice sites
///
typedef union
{
	struct
	{
		current_t evn;		//!< even lattice sites
		current_t odd;		//!< odd lattice sites
	};
	current_t evenodd[2];	//!< even and odd lattice sites stored in array
}
currsum_even_odd_t;


//_______________________________________________________________________________________________________________________
//

currsum_even_odd_t Field_Synchronize(const double time, field_t *f);


//_______________________________________________________________________________________________________________________
//

double Field_Elongation(const field_t *f);

double Field_Momentum(const field_t *f);

double Field_Energy(const field_t *f);



#endif
