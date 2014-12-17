/// \file current.h
/// \brief current J data structures, both physical and normal mode coordinates.
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

#ifndef CURRENT_H
#define CURRENT_H


//_______________________________________________________________________________________________________________________
///
/// \brief current J in physical coordinates
///
typedef struct
{
	double r;			//!< elongation current
	double p;			//!< momentum current
	double e;			//!< energy current
}
current_t;


//_______________________________________________________________________________________________________________________
///
/// \brief current J in normal modes
///
typedef struct
{
	double m[3];		//!< modes -1, 0 and 1
}
currnrm_t;



#endif
