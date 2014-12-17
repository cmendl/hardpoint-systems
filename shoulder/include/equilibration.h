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
	double hbeta;		//!< V_POT_HEIGHT * beta
	double Zspat;		//!< spatial part of the partition function
	double tcrit;		//!< position of discontinuity of derivative
	double coef1;		//!< coefficient 1 used in 'ElongationInverseCDF'
	double coef2;		//!< coefficient 2 used in 'ElongationInverseCDF'
	double sigma;		//!< 1/sqrt(beta), variance of the Gaussian distribution for the momentum
}
ensemble_params_t;

void EnsembleParams_Fill(const double pressure, const double beta, ensemble_params_t *params);


double AverageElongation(const double pressure, const double beta);

double AverageEnergy(const double pressure, const double beta);


void EquilibrateCanonical(const ensemble_params_t *params, randseed_t *seed, field_t *f);


//_______________________________________________________________________________________________________________________
//


void EquilibrateMicrocanonical(const double elongation, const double energy, randseed_t *seed, field_t *f);



#endif
