/// \file fft_transform.c
/// \brief Encapsulate FFTW plans and data arrays for input and output.
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

#include "fft_transform.h"
#include <memory.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Allocate memory for FFT transformation
///
void FFTtransform_Allocate(const unsigned int n, const int sign, fft_transform_t *trans)
{
	trans->fft_inp = (fftw_complex *)fftw_malloc(n*sizeof(fftw_complex));
	trans->fft_out = (fftw_complex *)fftw_malloc(n*sizeof(fftw_complex));

	trans->plan = fftw_plan_dft_1d(n, trans->fft_inp, trans->fft_out, sign, FFTW_ESTIMATE);

	// init with zeros
	memset(trans->fft_inp, 0, n*sizeof(fftw_complex));
	memset(trans->fft_out, 0, n*sizeof(fftw_complex));
}


//_______________________________________________________________________________________________________________________
///
/// \brief Delete FFT transformation structure (free memory)
///
void FFTtransform_Delete(fft_transform_t *trans)
{
	fftw_destroy_plan(trans->plan);

	fftw_free(trans->fft_out);
	fftw_free(trans->fft_inp);
}
