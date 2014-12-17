//	Equilibration of the lattice field variables: initialize elongations
//	and momenta randomly according to (micro-)canonical ensemble statistics.
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

#ifndef EQUILIBRATION_H
#define EQUILIBRATION_H

#include "field.h"
#include "random.h"


//_______________________________________________________________________________________________________________________
///
/// \brief Canonical ensemble parameters used for the fast evaluation of the
/// inverse cumulative distribution function (CDF) of the spatial part,
/// and the Gaussian distribution for the momentum
///
typedef struct
{
	double inv_q;		//!< 1 / (pressure * beta)
	double sigma[2];	//!< 1 / sqrt(mass_i * beta), variance of the Gaussian distribution for the velocity
}
ensemble_params_t;

void EnsembleParams_Fill(const double pressure, const double beta, const double mass[2], ensemble_params_t *params);

double ElongationInverseCDF(const ensemble_params_t *params, const double t);


double AverageElongation(const double pressure, const double beta);

double AverageEnergy(const double beta);


void EquilibrateCanonical(const double pressure, const double beta, randseed_t *seed, field_t *f);


//_______________________________________________________________________________________________________________________
//

double AveragePressure(const double elongation, const double energy);

double AverageBeta(const double energy);


void EquilibrateMicrocanonical(const double elongation, const double energy, randseed_t *seed, field_t *f);



#endif
