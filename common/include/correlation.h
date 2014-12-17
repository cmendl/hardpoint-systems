/// \file correlation.h
/// \brief Lattice point data and covariance matrix structure for calculating correlations in space and time.
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

#ifndef CORRELATION_H
#define CORRELATION_H


//_______________________________________________________________________________________________________________________
///
/// \brief Lattice point data (elongation, momentum and energy)
///
typedef union
{
	double g[3];		//!< vector of entries
	struct
	{
		double r;		//!< elongation
		double p;		//!< momentum
		double e;		//!< energy
	};
}
pointdata_t;


//_______________________________________________________________________________________________________________________
///
/// \brief Covariance matrix entries
///
typedef union
{
	// column-major order
	double mat[9];
	struct
	{
		double rr, pr, er;	// first column
		double rp, pp, ep;	// second
		double re, pe, ee;	// third
	} entries;
}
cov_mat_t;



#endif
