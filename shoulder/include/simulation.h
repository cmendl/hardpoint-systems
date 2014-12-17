//	"Sweep" the event table, i.e., process all active events within current time period.
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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "field.h"
#include "event_table.h"

//_______________________________________________________________________________________________________________________
//
// both the field and event table are modified

void CollectEvents(field_t *f, event_table_t *tab);

current_t SweepTable(field_t *f, event_table_t *tab);

void SweepTableJSite(field_t *f, event_table_t *tab, current_t J[NUM_SITES]);

void GlobalTimeshift(const double t_shift, field_t *f, event_table_t *tab);



#endif
