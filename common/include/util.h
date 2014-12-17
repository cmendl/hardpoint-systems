//	Utility functions.
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

#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>


//_______________________________________________________________________________________________________________________
///
/// \brief square function x -> x^2
///
static inline double square(const double x)
{
	return x*x;
}


//_______________________________________________________________________________________________________________________
///
/// \brief maximum of two numbers
///
static inline double maxf(const double x, const double y)
{
	if (x >= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

//_______________________________________________________________________________________________________________________
///
/// \brief minimum of two numbers
///
static inline double minf(const double x, const double y)
{
	if (x <= y)
	{
		return x;
	}
	else
	{
		return y;
	}
}


//_______________________________________________________________________________________________________________________
//


int ReadData(const char *filename, void *data, const size_t size, const size_t n);

int WriteData(const char *filename, const void *data, const size_t size, const size_t n, const bool append);



#endif
