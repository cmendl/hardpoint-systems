//	Random number generator based on George Marsaglia's Xorshift RNG algorithm improved by Sebastiano Vigna's xorshift1024*.
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

#ifndef RANDOM_H
#define RANDOM_H

#include <stdint.h>


//_______________________________________________________________________________________________________________________
///
/// \brief Random number generator based on George Marsaglia's Xorshift RNG algorithm
/// and improved by Sebastiano Vigna
///
typedef struct
{
	uint64_t s[16];		//!< seed consists of 16*64 = 1024 bits
	unsigned int p;		//!< counter
}
randseed_t;

void Random_SeedInit(uint64_t n, randseed_t *seed);


//_______________________________________________________________________________________________________________________
//

uint64_t Random_GetUint(randseed_t *seed);

double Random_GetUniform(randseed_t *seed);

double Random_GetNormal(randseed_t *seed);



#endif
