//	Calculate time correlations using FFT.
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

#ifndef TIME_CORRELATION_H
#define TIME_CORRELATION_H

#include "correlation.h"
#include "fft_transform.h"


//_______________________________________________________________________________________________________________________
///
/// \brief Time correlation structure
///
typedef struct
{
	pointdata_t *avr;			//!< averages
	cov_mat_t *sqr;				//!< averages of squares for each time shift
	unsigned int num_t;			//!< number of time points

	// use FFT for calculating correlations
	fft_transform_t fft_trans[2][3];	//!< forward FFT for fast convolution
	fft_transform_t fft_back;			//!< inverse FFT
}
time_correlation_t;


void TimeCorrelation_Allocate(const unsigned int num_t, time_correlation_t *tcorr);

void TimeCorrelation_Delete(time_correlation_t *tcorr);

void TimeCorrelation_Reset(time_correlation_t *tcorr);


//_______________________________________________________________________________________________________________________
//


void TimeCorrelation_Collect(const pointdata_t *data, time_correlation_t *tcorr);


void TimeCorrelation_ScalarMultiply(const double x, time_correlation_t *tcorr);



#endif
