//	Calculate spatial (lattice site) correlations, reference implementation for testing.
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

#ifndef SPATIAL_CORRELATION_REF_H
#define SPATIAL_CORRELATION_REF_H

#include "correlation.h"		// for 'pointdata_t'


//_______________________________________________________________________________________________________________________
///
/// \brief Spatial correlation structure
///
typedef struct
{
	pointdata_t avr[2];			//!< averages of elongation, momentum and energy at time 0 and t
	cov_mat_t *sqr;				//!< averages of squares, i.e., < r_0(0) r_j(t) >, < p_0(0) r_j(t) > ...
}
spatial_correlation_ref_t;


void SpatialCorrelationRef_Allocate(spatial_correlation_ref_t *scorr);

void SpatialCorrelationRef_Delete(spatial_correlation_ref_t *scorr);

void SpatialCorrelationRef_Reset(spatial_correlation_ref_t *scorr);


//_______________________________________________________________________________________________________________________
//


void SpatialCorrelationRef_Collect(const pointdata_t *data0, const pointdata_t *data1, spatial_correlation_ref_t *scorr);


void SpatialCorrelationRef_ScalarMultiply(const double x, spatial_correlation_ref_t *scorr);



#endif
