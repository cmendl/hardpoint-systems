/// \file random.c
/// \brief Random number generator based on George Marsaglia's Xorshift RNG algorithm improved by Sebastiano Vigna's xorshift1024*.
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

#include "random.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Initialize seeds by applying Sebastiano Vigna's xorshift64* to the number 'n'
///
void Random_SeedInit(uint64_t n, randseed_t *seed)
{
	// must be non-zero
	assert(n > 0);

	unsigned int i;
	for (i = 0; i < 16; i++)
	{
		n ^= n >> 12;	// a
		n ^= n << 25;	// b
		n ^= n >> 27;	// c

		// 2685821657736338717 = 72821711 * 36882155347, from Pierre L'Ecuyer's paper
		seed->s[i] = n * 2685821657736338717LL;
	}

	seed->p = 0;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Use George Marsaglia's Xorshift RNG algorithm
/// improved by Sebastiano Vigna's xorshift1024* using a multiplication of the output
/// to draw an unsigned integer; generator period is (2^1024 - 1)
///
/// References:
///   - George Marsaglia. Xorshift RNGs. Journal of Statistical Software 8, 1-6 (2003)
///   - Sebastiano Vigna. An experimental exploration of Marsaglia's xorshift generators, scrambled. (2014)
///   - Pierre L'Ecuyer. Tables of linear congruential generators of different sizes and good lattice structure. Mathematics of Computation 68, 249-260 (1999) 
///
uint64_t Random_GetUint(randseed_t *seed)
{
	uint64_t s0 = seed->s[seed->p];
	seed->p = (seed->p + 1) & 15;
	uint64_t s1 = seed->s[seed->p];
	s1 ^= s1 << 31;		// a
	s1 ^= s1 >> 11;		// b
	s0 ^= s0 >> 30;		// c
	seed->s[seed->p] = s0 ^ s1;

	// 1181783497276652981 = 769 * 13611541 * 112902689, from Pierre L'Ecuyer's paper
	return seed->s[seed->p] * 1181783497276652981LL;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Draw a uniform random sample from the open interval (0, 1)
///
double Random_GetUniform(randseed_t *seed)
{
	// 0 <= u < 2^64
	uint64_t u = Random_GetUint(seed);
	// factor is 1/(2^64 + 2), such that result is strictly between 0 and 1
	return (u + 1.0) * 5.4210108624275221694495168e-20;
}


//_______________________________________________________________________________________________________________________
///
/// \brief Draw a normal (Gaussian) random sample with mean 0 and standard deviation 1
///
double Random_GetNormal(randseed_t *seed)
{
	// use Box-Muller transform
	double u1 = Random_GetUniform(seed);
	double u2 = Random_GetUniform(seed);
	double r = sqrt(-2.0*log(u1));
	return r * sin(2.0*M_PI * u2);
}
